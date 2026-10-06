from pathlib import Path


class PathManager:
    def __init__(self):
        this_file = Path(__file__).resolve()
        root_project = this_file.parents[2]
        root = Path(
            '/exchange/healthds/pQTL/BELIEVE')
        if root_project.parent != root:
            print(f"Project is outside {root}; using {root_project}")
            print("Using root project: {}".format(root_project))
            root = root_project

        literature_table_path = 'literature_table'
        literature_config_path = 'literature_config'
        literature_files_path = 'literature_files'
        literature_harmonized_path = 'literature_harmonized'
        literature_gwasstudio_files_path = 'literature_gwasstudio_files'
        literature_gwasstudio_output_path = 'literature_gwasstudio_output'

        self.inputs = {
            'literature_table_raw' : Path(root_project, literature_table_path, 'literature_table_all_somalogic.xlsx'),
            'literature_table' : Path(root_project, literature_table_path, 'literature_table_all_somalogic_allstudies.xlsx'),
            'literature_table_cleaned' : Path(root_project, literature_table_path, 'literature_table_all_somalogic_cleaned.xlsx'),
            'literature_table_harmonized' : Path(root_project, literature_table_path, 'literature_table_all_somalogic_harmonized.xlsx'),
        }
        self.config = {
            'config_harmonize_build38' : Path(root_project, literature_config_path, 'config_harmonize_build38.yml'),
            'config_harmonize_build37' : Path(root_project, literature_config_path, 'config_harmonize_build37.yml'),
            'config_harmonize_build37_bcftools' : Path(root_project, literature_config_path, 'liftover_test/config_harmonize_build37_bcftools.yml'),
            'config_harmonize_build37_gwaslab_bridge' : Path(root_project, literature_config_path, 'liftover_test/config_harmonize_build37_gwaslab_bridge.yml'),
            'config_harmonize_build37_gwaslab_standard' : Path(root_project, literature_config_path, 'liftover_test/config_harmonize_build37_gwaslab_standard.yml'),
            'believe_metadata' : Path(root_project, literature_config_path, 'believe_metadata.tsv'),
            'literature_panel' : Path(root_project, literature_config_path, 'literature_protein_panel.tsv'),
            'panels_map' : Path(root_project, literature_config_path, 'believe_literature_panels_map.tsv'),
        }
        self.files = {
            'pqtl_sun_ukb_csa' : Path(root_project, literature_files_path, 'sun_ukb_st11.csv'),
            'pqtl_interval_chris_meta' : Path(root_project, literature_files_path, 'interval_chris_meta_st3.csv'),
            'pqtl_decode_2023' : Path(root_project, literature_files_path, 'pqtl_decode_2023/pqtl_decode_2023_leadsnps.csv'),
            'pqtl_CKB_SomaScan' : Path(root_project, literature_files_path, 'ckb_somascan.csv'),
            'pqtl_CKB_Olink' : Path(root_project, literature_files_path, 'ckb_olink.csv'),
        }
        self.outputs = {
            'literature_harmonized': Path(root, literature_harmonized_path),
            'literature_gwasstudio_files': Path(root, literature_gwasstudio_files_path),
            'literature_gwasstudio_output': Path(root, literature_gwasstudio_output_path),
            
            # Liftover test outputs
            'literature_table_harmonized_bcftools' : Path(root, literature_harmonized_path, 'liftover_test/bcftools/literature_table_all_somalogic_harmonized.xlsx'),
            'literature_table_liftover_bcftools' : Path(root, literature_harmonized_path, 'liftover_test/bcftools/literature_table_all_somalogic_liftover.xlsx'),
            'literature_table_harmonized_gwaslab_bridge' : Path(root, literature_harmonized_path, 'liftover_test/gwaslab_bridge/literature_table_all_somalogic_harmonized.xlsx'),
            'literature_table_harmonized_gwaslab_standard' : Path(root, literature_harmonized_path, 'liftover_test/gwaslab_standard/literature_table_all_somalogic_harmonized.xlsx'),
        }

    def get_inputs(self):
        return self.inputs

    def get_config(self):
        return self.config

    def get_files(self):
        return self.files

    def get_outputs(self):
        return self.outputs

    def get_output(self, label, exists=True):
        output_path = self.outputs.get(label, None)
        if exists and not self.outputs[label].exists():
            return None
        return output_path

