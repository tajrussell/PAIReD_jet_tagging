import os
import math
import torch
from torch import Tensor
from weaver.utils.logger import _logger
from weaver.utils.import_tools import import_module

ParticleTransformerTagger = import_module(
    os.path.join(os.path.dirname(__file__), 'ParticleTransformer.py'), 'ParT').ParticleTransformerTagger


def get_model(data_config, **kwargs):

    cfg = dict(
        pf_input_dim=len(data_config.input_dicts['pf_features']),
        sv_input_dim=len(data_config.input_dicts['sv_features']),
        # network configurations
        pair_input_dim=4,
        num_classes=3*len(data_config.label_value),
        embed_dims=[128, 512, 128],
        pair_embed_dims=[64, 64, 64],
        num_heads=8,
        num_layers=8,
        num_cls_layers=2,
        block_params=None,
        cls_block_params={'dropout': 0, 'attn_dropout': 0, 'activation_dropout': 0},
        fc_params=[],
        activation='gelu',
        # misc
        trim=True,
        for_inference=False,
    )
    #kwargs.pop('loss_gamma')
    cfg.update(**kwargs)
    _logger.info('Model config: %s' % str(cfg))

    model = ParticleTransformerTagger(**cfg)

    model_info = {
        'input_names': list(data_config.input_names),
        'input_shapes': {k: ((1,) + s[1:]) for k, s in data_config.input_shapes.items()},
        'output_names': ['softmax'],
        'dynamic_axes': {**{k: {0: 'N', 2: 'n_' + k.split('_')[0]} for k in data_config.input_names}, **{'softmax': {0: 'N'}}},
    }

    return model, model_info


class LogCoshLoss(torch.nn.L1Loss):
    __constants__ = ['reduction']

    def __init__(self, reduction: str = 'mean') -> None:
        super(LogCoshLoss, self).__init__(None, None, reduction)

    def forward(self, input: Tensor, target: Tensor) -> Tensor:
        #print('LogCoshLoss steps')
        x = (input - target)
        loss = (x + torch.nn.functional.softplus(-2. * x) - math.log(2))
        #print('loss', loss)
        if self.reduction == 'none':
            return loss
        elif self.reduction == 'mean':
            return loss.mean()
        elif self.reduction == 'sum':
            return loss.sum()
        

class QuantileLoss(torch.nn.L1Loss):
    __constants__ = ['reduction']

    def __init__(self, reduction: str = 'mean', quantile=0.16) -> None:
        super(QuantileLoss, self).__init__(None, None, reduction)
        self.quantile = quantile

    def forward(self, input: Tensor, target: Tensor) -> Tensor:
        #print('Quantile Loss steps')
        z = (target - input)
        #print('z', z)
        loss = (self.quantile * z * (z>=0) + (self.quantile - 1) * z * (z<0))
        #print('loss', loss)
        if self.reduction == 'none':
            return loss
        elif self.reduction == 'mean':
            return loss.mean()
        elif self.reduction == 'sum':
            return loss.sum()


class RegressionLoss(torch.nn.L1Loss):
    __constants__ = ['reduction']

    def __init__(self, reduction: str = 'mean', factor_reg=1., factor_err=1.) -> None:
        super(RegressionLoss, self).__init__(None, None, reduction)
        self.loss_reg_fn = LogCoshLoss()
        self.loss_err_fn_minus = QuantileLoss(quantile=0.16)
        self.loss_err_fn_plus = QuantileLoss(quantile=0.84)
        self.factor_reg = factor_reg
        self.factor_err = factor_err

    def forward(self, preds: Tensor, target_reg: Tensor) -> Tensor:
        #print(target_reg.shape)
        #print(preds.shape)
        input_reg = preds[:, 0]
        input_err_plus = preds[:, 1]
        input_err_minus = preds[:, 2]
        #print(input_reg.shape)
        loss_reg = self.loss_reg_fn(input_reg, target_reg)
        loss_err = self.loss_err_fn_plus(input_err_plus, target_reg) + self.loss_err_fn_minus(input_err_minus, target_reg)
        loss = self.factor_reg * loss_reg + self.factor_err * loss_err
        return loss, {'reg': loss_reg.item(), 'err': loss_err.item()}


def get_loss(data_config, **kwargs):
    factor_reg = 1.0#kwargs.get('factor_reg', 1)
    factor_err = 1.0#kwargs.get('factor_err', 1)
    return RegressionLoss(factor_reg=factor_reg, factor_err=factor_err)

def get_metrics():
    def mse_median(y_true, y_pred):
        return mean_squared_error(y_true, y_pred[:, 1])

    def mae_median(y_true, y_pred):
        return mean_absolute_error(y_true, y_pred[:, 1])

    def coverage(y_true, y_pred):
        y_lower, y_median, y_upper = y_pred.T
        return np.mean((y_true >= y_lower) & (y_true <= y_upper))

    def interval_width(y_true, y_pred):
        y_lower, _, y_upper = y_pred.T
        return np.mean(y_upper - y_lower)

    return {
        "mse_median": mse_median,
        "mae_median": mae_median,
        "coverage": coverage,
        "interval_width": interval_width,
    }
