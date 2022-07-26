#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


char* make_def_path(const char* name)
{
    int pid;
    char procpath[25], exepath[1024]; // 1024 SHOULD be enough (ha!)
    char *dirpath, *filepath, *exedir;
    ssize_t linklen, destlen;

    pid = getpid();
    sprintf(procpath, "/proc/%d/exe", pid);
    linklen = readlink(procpath, exepath, 1024);
    if(linklen < 0) {
        perror("lstat");
        exit(EXIT_FAILURE);
    }

    if(linklen > 1023) {
        fprintf(stderr, "Executable path changed or too long. This should not occur.\n");
        exit(EXIT_FAILURE);
    }

    exepath[linklen + 1] = '\0';
    exedir = strdup(exepath);
    dirpath = dirname(exedir);
    destlen = strlen(dirpath) + 1 + strlen(name);
    filepath = malloc(destlen + 1);
    if(filepath == NULL) {
        fprintf(stderr, "Insufficient memory to construct path\n");
        exit(EXIT_FAILURE);
    }

    sprintf(filepath, "%s/%s", dirpath, name);

    free(dirpath);
    return filepath;
}

