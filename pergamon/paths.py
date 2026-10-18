"""Repository-local runtime paths for Pergamon."""

import os

import tdpy
from tdpy.paths import RepositoryPaths


PATH_ENV_VAR = "PERGAMON_PATH"
_REPOSITORY_PATHS = RepositoryPaths(PATH_ENV_VAR)

get_repository_path = _REPOSITORY_PATHS.get_repository_path
get_data_path = _REPOSITORY_PATHS.get_data_path
get_visuals_path = _REPOSITORY_PATHS.get_visuals_path


def retr_pathgaia(pathbase=None):
	"""Return Pergamon base, data, and image paths for Gaia workflows."""

	if pathbase is None:
		pathbase = tdpy.ensr_path(get_repository_path())
	else:
		pathbase = tdpy.ensr_path(pathbase)

	return {
		'pathbase': pathbase,
		'pathdata': tdpy.ensr_path(os.path.join(pathbase, 'data')),
		'pathimag': tdpy.ensr_path(os.path.join(pathbase, 'visuals')),
	}