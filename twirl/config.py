from pathlib import Path
import yaml
from yaml import Loader
import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table

info = print
package_name = "twirl"


class ConfigManager:
    def __init__(self):

        self.config = None

        self.folder_path = Path.home() / f".{package_name}"
        self.folder_path.mkdir(exist_ok=True)

        self.config_file = self.folder_path / "config"

        self.check_config_file(load=True)

        self.catalog_path = self.folder_path / f"gaia_14.fits.gz"
        self.current_catalog = None

    def check_config_file(self, load=False):

        if self.config_file.exists():
            with self.config_file.open(mode="r") as file:
                if load:
                    self.config = yaml.load(file.read(), Loader=Loader)
        else:
            info(f"A config file as been created in {self.folder_path}")
            self.config = {"color": "blue"}
            with self.config_file.open(mode="w") as file:
                yaml.dump(self.config, file, default_flow_style=False)

    def save(self):
        self.check_config_file()
        with self.config_file.open(mode="w") as file:
            yaml.dump(self.config, file, default_flow_style=False)

    def download_catalog(self):
        info("downloading lt14 catalog (~300Mb)")
        catalog = requests.get("https://github.com/lgrcia/twirl-catalogs/raw/master/lt14gmag-result.fits.gz").content
        self.catalog_path.open(mode="wb").write(catalog)  
        self.load_catalog()

    def load_catalog(self):
        if not self.catalog_path.exists():
            self.download_catalog()
        
        if self.current_catalog is None:
            info("loading catalog")
            gaia_sources = Table.read(self.catalog_path, hdu=1)
            self.current_catalog = np.array([gaia_sources["ra"].data, gaia_sources["dec"].data]).T

