def get_config(cfg, key, default=None):
    if isinstance(cfg, dict):
        if key in cfg and cfg[key] is not None:
            return cfg[key]
        for v in cfg.values():
            result = get_config(v, key, default)
            if result is not None:
                return result
    elif isinstance(cfg, list):
        for item in cfg:
            result = get_config(item, key, default)
            if result is not None:
                return result
    return default