import uproot
import awkward as ak

def read_weights_from_file(file_path, reweight_classes):
    weights = dict()
    with uproot.open(file_path) as file:
        tree = file['tree']
        for label in reweight_classes:
            weights[label] = tree['weights_'+label].array(library='ak')
    return weights

def write_weights_to_text(weights, file_path, dijet_pt_bins, mc_mass_bins, reweight_classes, class_weights, reweight_threshold):
    with open(file_path, 'w') as file:
        file.write('weights:\n')
        file.write('  use_precomputed_weights: false\n')
        file.write('  reweight_method: ref\n')
        file.write('  reweight_vars:\n')
        file.write('    dijet_pt:\n')
        for pt in dijet_pt_bins:
            file.write(f'    - {pt}\n')
        file.write('    MC_gendijet_mass:\n')
        for mass in mc_mass_bins:
            file.write(f'    - {mass}\n')
        file.write(f'  reweight_threshold: {reweight_threshold}\n')
        file.write('  reweight_classes:\n')
        for label in reweight_classes:
            file.write(f'  - {label}\n')
        file.write('  class_weights:\n')
        for weight in class_weights:
            file.write(f'  - {weight}\n')
        file.write('  reweight_hists:\n')

        for label in reweight_classes:
            file.write(f'    {label}:\n')
            label_weights = weights[label]
            for pt_idx in range(len(dijet_pt_bins) - 1):
                file.write('    -')
                for mass_idx in range(len(mc_mass_bins) - 1):
                    weight_value = label_weights[pt_idx*(len(mc_mass_bins)-1) + mass_idx]
                    if mass_idx == 0:
                        file.write(f' - {weight_value}\n')
                    else:
                        file.write(f'      - {weight_value}\n')

if __name__ == '__main__':
    reweight_classes = ['label_BB', 'label_CC', 'label_bb', 'label_bx', 'label_cx', 'label_ll']
    weights = read_weights_from_file('reweight.root', reweight_classes)
    dijet_pt_bins = [0.00000000e+00, 3.03419350e+02, 4.66926422e+02, 9.08118271e+02, 9.99999900e+06]
    mc_mass_bins = [25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245]
    class_weights = [1.0, 1.0, 0.12, 0.22, 0.85]
    reweight_threshold = 0.1
    write_weights_to_text(weights, 'weights.txt', dijet_pt_bins, mc_mass_bins, reweight_classes, class_weights, reweight_threshold)
    print('Weights saved to weights.txt')