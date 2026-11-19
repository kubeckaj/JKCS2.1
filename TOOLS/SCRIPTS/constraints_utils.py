import os

TRAIN_CONF_FILE = 'gen_config.yaml'
PARENT_FOLDER = '../../'

def real_path(path):
    if path == None or path == '':
        return False
    
    if os.path.exists(path):
        return path
    
    path = PARENT_FOLDER+path
    if os.path.exists(path):
        return path
    
    return False

def save_config_file(args, ckpt_folder):
    import yaml

    KEYS_TO_SAVE = [
        "config_root",
        "model",
        "use_molgraphconvfeat",
        "enc_monomers",
        "dataset",
        "data_path",
        "max_flexibility",
        "include_mono",
        "keep_h",
    ]

    cfg = {k: getattr(args, k) for k in KEYS_TO_SAVE}
    cfg["load_ckpt"] = os.path.basename(ckpt_folder)
    with open(os.path.join(ckpt_folder, TRAIN_CONF_FILE), 'w') as f:
        yaml.safe_dump(cfg, f)

    return

def load_config_file(cfg):
    import yaml
    
    ckpt_folder = cfg["model_ckpt_folder"]
    with open(os.path.join(ckpt_folder, TRAIN_CONF_FILE), 'r') as f:
        cfg_ckpt = yaml.safe_load(f)

    cfg_ckpt["load_ckpt"] = os.path.basename(ckpt_folder)
    cfg.update(cfg_ckpt)
    return cfg
