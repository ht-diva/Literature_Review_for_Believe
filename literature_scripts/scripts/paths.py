from pathlib import Path


class PathManager:
    def __init__(self):
        root = Path(
            '/exchange/healthds/pQTL/BELIEVE')
        if not root.exists():
            exit("Path not found: {}".format(root))
        else:
            self._root_path = root
        this_file = Path(__file__).resolve()
        root_project = this_file.parents[2]

        literature_table_path = 'literature_table'
        literature_config_path = 'literature_config'
        literature_files_path = 'literature_files'
        literature_harmonized_path = 'literature_harmonized'
        literature_gwasstudio_files_path = 'literature_gwasstudio_files'
        gwasstudio_output_path = 'gwasstudio_output'

        self.inputs = {
            'literature_table_raw' : Path(root, literature_table_path, 'literature_table_all_somalogic.xlsx'),
            'literature_table' : Path(root, literature_table_path, 'literature_table_all_somalogic_cleaned.xlsx'),
        }
        self.config = {
            'config_harmonize' : Path(root_project, literature_config_path, 'config_harmonize.yml'),
            'believe_metadata' : Path(root_project, literature_config_path, 'believe_metadata.tsv'),
        }
        self.files = {
            'pqtl_sun_ukb_csa' : Path(root_project, literature_files_path, 'sun_ukb_st11.csv'),
            'pqtl_interval_chris_meta' : Path(root_project, literature_files_path, 'interval_chris_meta_st3.csv'),
            'pqtl_decode_2023' : Path(root_project, literature_files_path, 'pqtl_decode_2023/pqtl_decode_2023_leadsnps.csv'),
        }
        self.outputs = {
            'literature_harmonized': Path(root, literature_harmonized_path),
            'literature_gwasstudio_files': Path(root, literature_gwasstudio_files_path),
            'gwasstudio_output': Path(root, gwasstudio_output_path),
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
