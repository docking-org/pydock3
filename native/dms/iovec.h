/* struct iovec from <sys/uio.h>, and caddr_t (pydock3: which MinGW does not have) */
#ifndef IOVEC_H
#define IOVEC_H
#ifdef _WIN32
#include <stddef.h>
typedef char	*caddr_t;
struct iovec {
	void	*iov_base;
	size_t	iov_len;
};
#else
#include <sys/types.h>
#include <sys/uio.h>
#endif
#endif
