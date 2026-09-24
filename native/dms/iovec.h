/* struct iovec, from <sys/uio.h> (pydock3: which MinGW does not have) */
#ifndef IOVEC_H
#define IOVEC_H
#ifdef _WIN32
#include <stddef.h>
struct iovec {
	void	*iov_base;
	size_t	iov_len;
};
#else
#include <sys/uio.h>
#endif
#endif
