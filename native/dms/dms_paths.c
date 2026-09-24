/*
 * Paths of files installed next to the dms executable, e.g. the dmsd server.
 * (pydock3: found from argv[0] instead of /proc/<pid>/exe, which only Linux has;
 * pydock3 runs dms by its full path.)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char	*exe_dir = ".";

void set_exe_dir(const char *argv0)
{
	const char	*sep = strrchr(argv0, '/');
#ifdef _WIN32
	const char	*backslash = strrchr(argv0, '\\');

	if (backslash != NULL && (sep == NULL || backslash > sep))
		sep = backslash;
#endif
	if (sep != NULL) {
		exe_dir = malloc(sep - argv0 + 1);
		if (exe_dir == NULL) {
			fprintf(stderr, "Insufficient memory to construct path\n");
			exit(EXIT_FAILURE);
		}
		memcpy(exe_dir, argv0, sep - argv0);
		exe_dir[sep - argv0] = '\0';
	}
}

char *make_def_path(const char *name)
{
	char	*path = malloc(strlen(exe_dir) + 1 + strlen(name) + 1);

	if (path == NULL) {
		fprintf(stderr, "Insufficient memory to construct path\n");
		exit(EXIT_FAILURE);
	}
	sprintf(path, "%s/%s", exe_dir, name);
	return path;
}
