import json, os
import logging
from pathlib import Path

def load_config(CONFIG_FILE):
    # Get the absolute path of the current script
    script_dir = Path(__file__).resolve().parent
    # Look for the config file in the parent directory of the script
    config_path = script_dir.parent / CONFIG_FILE
    logging.info(f"Looking for config file at: {config_path}")
    try:
        with open(config_path, 'r') as config_file:
            config = json.load(config_file)
        logging.info(f"Config file loaded successfully from {config_path}")
        return config
    except FileNotFoundError:
        logging.warning(f"Config file {CONFIG_FILE} not found at {config_path}. Using default values.")
        return {}
    except json.JSONDecodeError:
        logging.error(f"Config file {CONFIG_FILE} is not valid JSON. Using default values.")
        return {}
    
def initialize_environment(CONFIG_FILE):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Log the current working directory
    logging.info(f"Current working directory: {os.getcwd()}")
    
    config = load_config(CONFIG_FILE)
    local_env = config.get('localenv', '')
    
    env_args = {}
    
    if not local_env:
        logging.warning("Local environment path not set in config. Assuming system-wide installation.")
        # Check for samtools in system PATH
        samtools_path = 'samtools'
        env_args["samtools"] = samtools_path
        
        # Check for Python (for pysam)
        python_path = 'python'
        env_args["pysam"] = python_path
    else:
        # Check for samtools
        samtools_path = os.path.join(local_env, 'bin', 'samtools')
        if os.path.exists(samtools_path):
            env_args["samtools"] = samtools_path
        else:
            logging.warning("Samtools not found in the local environment. Please install samtools.")

        # Check for pysam (assuming it's installed in the Python environment)
        python_path = os.path.join(local_env, 'bin', 'python')
        if os.path.exists(python_path):
            # We'll use the Python interpreter path, as pysam is a Python package
            env_args["pysam"] = python_path
        else:
            logging.warning("Python not found in the local environment. Please ensure Python and pysam are installed.")

    if not env_args:
        logging.error("Neither samtools nor Python (for pysam) were found.")
        logging.error("Please install the required software or specify the correct local environment path in the config file.")
    else:
        logging.info("Environment initialized successfully.")
    
        env_vars_file = Path(__file__).parent.parent / 'utils' / 'env_vars.json'
        
        # Create the utils directory if it doesn't exist
        env_vars_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Write the environment variables to the JSON file
        with open(env_vars_file, 'w') as f:
            json.dump(env_args, f, indent=4)
        
        logging.info(f"Environment variables written to {env_vars_file}")

def add_parser(subparsers):
    init_parser = subparsers.add_parser('init', help='Initialize the Atlas environment')
    init_parser.add_argument('-c', '--config', type=str, default='config.json', help='Path to the config file')
    return init_parser
