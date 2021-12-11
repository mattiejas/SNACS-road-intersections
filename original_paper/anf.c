/***************************************************************************

                                A N F

                     Approximate Neighbourhood Function


   This source code is (c) Copyright 2001 by Christopher R. Palmer.
   It may not be redistributed without the express consent of the
   author.


   -------------------------------------------------------------------------

   Building this program
   =====================

   To build this program, you can run:

        gcc -O anf.c -o anf -lpthread -lm

   ------------------------------------------------------------------------

   Basic Input / Output
   ====================

   With no arguments, this program will reads stdin with the following
   file format.  The first line contains the number of nodes in the
   graph, N.  The following lines contain pairs of nodes which are
   each in the range 0..N-1:

     x y

   which indicates a directed edge from node x to node y.  For example:

   6
   2 3
   3 4
   4 5
   5 2

   is an input that corresponds to the graph:

           0    2 ----> 3
                ^       |
                |       |
                |       ^
           1    5 <---- 4

   with two disconnected nodes, 0 and 1.  The output is a list of
   pairs

     h N(h)

   which indicates that the neighbourhood function at h is approximately
   N(h).  For example, one run of the previous data set produced

     1 13.727006
     2 17.019614
     3 21.636301
     4 21.636301

   which is an approximation of the neighbourhood function
   (the true value for this example is: N(1)=10, N(2)=18, N(>=3)=22).

   --------------------------------------------------------------------

   Command Line Arguments
   ======================

   This list includes only the supported command line arguments.  There
   are other arguments that are accepted which are not officially
   supported.  You may guess their meaning at your own risk...

   The following are general arguments:

     -mem n   : use up-to n MB of memory for the processing [default = 64MB]
     -verbose : print verbose information about the steps of the computation
     -n n     : the number of nodes, SEE BELOW!

   The following arguments modify the parameters of the algorithm:

     -k k     : average over k parallel approximations (larger k involves
                more time and memory but better accuracy) [default = 32]
     -bits b  : use log(n) + b bits for the estimates.  Larger values
                add variance, smaller values add truncation errors. [default=5]
     -it it   : only run up to it iterations, negative => infinity [default=-1]
     -start set : compute only for paths that start in this set
     -probe set : and go to nodes in this set.  if you specify probe with
		  no start set then start = probe (and visa-versa).

   The following arguments control the files used by the program:

     -past fname : use fname as the file to store the data for the past
                   iteration.  use this if your /tmp partition does not
                   contain enough space for all the data.  [default = n/a]
     -cur fname  : same, but for the current iteration     [default = n/a]
     -compress   : do leading ones compression on the bitmasks

   The -n argument is very strange.  When you specify the number of nodes
   on the command line, the input format is changed.  The file will then
   contain a pair of nodes on each input line (the first input line is
   not read).

   --------------------------------------------------------------------

   Version: 0.4
   Date:    October 10, 2001

   (c) Copyright 2001, Christopher R. Palmer
       crpalmer@cs.cmu.edu

 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>

#if EXPERIMENTAL_DISTANCE
#include "anf-dist.h"
#endif

/* do the FM stochastic estimation */
#define STOCHASTIC 0
#define SINGLE_THREADED 0

#if !SINGLE_THREADED
#include <pthread.h>
#endif

#define T ((unsigned long)time(NULL))

#define ULB(a) (((unsigned long)(a)) & (sizeof(unsigned long) - 1))

#define wordcopy(_dst, _src, _nbytes)                      \
	do                                                     \
	{                                                      \
		unsigned char *dst = (unsigned char *)_dst;        \
		unsigned char *src = (unsigned char *)_src;        \
		int nbytes = _nbytes;                              \
		int pre = sizeof(unsigned long) - ULB(dst);        \
		int i;                                             \
                                                           \
		if (ULB(src) != ULB(dst))                          \
		{                                                  \
			memcpy(dst, src, nbytes);                      \
			break;                                         \
		}                                                  \
                                                           \
		if (pre == sizeof(unsigned long))                  \
			pre = 0;                                       \
		if (pre > nbytes)                                  \
			pre = nbytes;                                  \
		for (i = 0; i < pre; i++)                          \
			*dst++ = *src++;                               \
		nbytes -= pre;                                     \
		while (nbytes >= sizeof(unsigned long))            \
		{                                                  \
			*(unsigned long *)dst = *(unsigned long *)src; \
			dst += sizeof(unsigned long);                  \
			src += sizeof(unsigned long);                  \
			nbytes -= sizeof(unsigned long);               \
		}                                                  \
		while (nbytes)                                     \
		{                                                  \
			*dst++ = *src++;                               \
			nbytes--;                                      \
		}                                                  \
	} while (0)

int verbose = 0;

#if !SINGLE_THREADED
typedef struct
{
	pthread_mutex_t m;
	pthread_cond_t c;
	int x;
} sem_t;
#endif

typedef struct
{
#if !SINGLE_THREADED
	pthread_t th;
	sem_t start;
	sem_t reply;
#endif
	enum
	{
		START,
		WORK,
		REPLY
	} state;
	int fd;
	struct
	{
		int is_read;
		void *buf;
		int off;
		int len;
		int *clen;
	} * r;
	int nr;
	int ar;
	long res;
	long total_len;
} iothread_t;

#define N_IN(nbuf) \
	((n + (nbuf)-1) / (nbuf))

#define BUCKET(index, nbucks, n) \
	((index) / N_IN(nbucks))

#define BIAS (1 + .31 / k)

int BYTES1;

#define BYTES(nent) \
	(((k * (logn + bits) + select_bits + compression_bits + 7) / 8) * (nent))
#define BYTES_IN(nbuf) \
	(BYTES(N_IN(nbuf)))

#define FILE_OFFSET(index, nbucks, n) \
	(BUCKET(index, nbucks, n) * BYTES_IN(nbucks))

typedef struct
{
	int u, v;
} edge_t;

int mem = 64; /* MB of memory to use */
int bits = 5; /* extra bits per node */
int k = 32;	  /* number of trials for each of the splits */
int nk = 1;	  /* number of k splits */
int do_compression = 0;
int compression_bits = 0; /* number of bits used for leading 1 counter */
int select_bits = 0;
char *pastfname = NULL;
char *curfname = NULL;
char *edgefname = NULL;
int *past_clens, *cur_clens;
int edges_bytes_avail = 0;
#if EXPERIMENTAL_DISTANCE
anf_dist_t *dist = NULL;
#endif

int n; /* number of nodes */
int logn;

unsigned char *buffers[5];

FILE **edgef; /* table of edge files */
int pastfd = -1;
int curfd = -1;

int nedgef; /* number of files */
int iedgef; /* currently accessed file */

int n_buckets[2]; /* number of equal sized buckets in past / cur bufs */
int *pastind;
int *curind;

edge_t *edges = NULL;
int total_edges;

#if !SINGLE_THREADED

int sem_create(sem_t *sem, int value)
{
	if (pthread_mutex_init(&sem->m, NULL) == -1)
	{
		perror("sem_create(mutex_init)");
		exit(1);
	}

	if (pthread_cond_init(&sem->c, NULL) == -1)
	{
		perror("sem_create(cond_init)");
		exit(1);
	}

	sem->x = value;

	return 0;
}

int sem_destroy(sem_t *sem)
{
	pthread_cond_destroy(&sem->c);
	pthread_mutex_destroy(&sem->m);
	free(sem);

	return 0;
}

int sem_P(sem_t *sem)
{
	pthread_mutex_lock(&sem->m);
	while (sem->x <= 0)
		pthread_cond_wait(&sem->c, &sem->m);
	sem->x--;
	pthread_mutex_unlock(&sem->m);

	return 0;
}

int sem_V(sem_t *sem)
{
	pthread_mutex_lock(&sem->m);
	if (sem->x++ == 0)
		pthread_cond_signal(&sem->c);
	pthread_mutex_unlock(&sem->m);

	return 0;
}

#endif

static void
iothread_doit(iothread_t *me)
{
	int i, res;
	int len;

	me->res = 0;

	for (i = 0; i < me->nr; i++)
	{
		if (me->r[i].off >= 0)
		{
			if ((res = lseek(me->fd, me->r[i].off, SEEK_SET)) < 0)
			{
				me->res = -errno;
				break;
			}
		}

		if (compression_bits && me->r[i].clen)
			len = *me->r[i].clen;
		else
			len = me->r[i].len;

		if (me->r[i].is_read)
		{
			res = read(me->fd, me->r[i].buf + (me->r[i].len - len), len);
		}
		else
		{
			res = write(me->fd, me->r[i].buf, len);
		}

		if (res < 0)
		{
			me->res = -errno;
			break;
		}
		else if (res == 0)
			break;
		else
			me->res += res;
	}
}

#if !SINGLE_THREADED

void *
iothread_main(void *input)
{
	iothread_t *me = input;

	for (;;)
	{
		sem_P(&me->start);
		assert(me->state == WORK);
		iothread_doit(me);
		me->state = REPLY;
		sem_V(&me->reply);
	}
}

#endif

iothread_t iot[3];

void iothread_create(iothread_t *io)
{
	io->state = START;
	io->ar = 2;
	io->r = malloc(sizeof(*io->r) * io->ar);

#if !SINGLE_THREADED

	sem_create(&io->start, 0);
	sem_create(&io->reply, 0);

	if (pthread_create(&io->th, NULL, iothread_main, io) == -1)
	{
		perror("pthread_create");
		exit(1);
	}
#endif
}

static void
dump_iothread(iothread_t *io)
{
	int i;

	fprintf(stderr, "fd %d ", io->fd);
	for (i = 0; i < io->nr; i++)
	{
		fprintf(stderr, "%s%s %d/%d bytes at %d to/from %p", i ? ", " : "", io->r[i].is_read ? "read" : "write", io->r[i].clen ? *io->r[i].clen : 0, io->r[i].len, io->r[i].off, io->r[i].buf);
	}
	fprintf(stderr, "\n");
}

/* arg 2+4i: void * buffer[i]  : anf buffer start
       3+4i: bool   is_read[i] : operation: read/write
       4+4i: int    off[i]     : offset in the file
       5+4i: int    len[i]     : true length of the buffer
       6+5i: int  * clen[i]    : compressed length of the buffer
                                 read for a read request
                                 written for a write request
       7+5n: NULL

   For each buffer that will be written, if we have enabled leading
   ones compression then in place transform the buffer into a compressed
   representation.  The representation is:

       (#ones-1)(otherbits-1) ... (other bits-n)(#ones-n) (blank space)
         byte 1   bytes 2 - ...

   where we always preserve whole bytes.  That is, we do not compress
   down to the bit level, instead compressing only the bytes that contain
   only leading ones.

   After an I/O operation completes, the wait function should uncompress
   all of the buffers (both read and write).

 */

void iothread_start(iothread_t *io, int fd, ...)
{
	va_list ap;
	void *buf;

	assert(io->state == START);

	io->total_len = 0;
	va_start(ap, fd);

	io->fd = fd;

	for (io->nr = 0; (buf = va_arg(ap, void *)) != NULL; io->nr++)
	{
		if (io->nr >= io->ar)
		{
			io->ar *= 2;
			io->r = realloc(io->r, sizeof(*io->r) * io->ar);
		}
		io->r[io->nr].buf = buf;
		io->r[io->nr].is_read = va_arg(ap, int);
		io->r[io->nr].off = va_arg(ap, int);
		io->r[io->nr].len = va_arg(ap, int);
		io->r[io->nr].clen = va_arg(ap, int *);

		if (!io->r[io->nr].is_read && compression_bits && io->r[io->nr].clen)
		{
			int i, j;
			char *bf = buf;

			for (i = j = 0; i < io->r[io->nr].len; i += BYTES1)
			{
				int compressed_bytes;

				bf[j++] = bf[i];
				compressed_bytes = ((bf[i] & ((1 << compression_bits) - 1)) * k - compression_bits - select_bits) / 8;
				if (!compressed_bytes)
					compressed_bytes = 1;
#if 1
				wordcopy(&bf[j], &bf[i + compressed_bytes], BYTES1 - compressed_bytes);
				j += BYTES1 - compressed_bytes;
#else
				for (l = compressed_bytes; l < BYTES1; l++)
					bf[j++] = bf[i + l];
#endif
			}
			*io->r[io->nr].clen = j;
			io->total_len += j;
		}
		else if (compression_bits && io->r[io->nr].clen)
		{
			io->total_len += *io->r[io->nr].clen;
		}
		else
		{
			io->total_len += io->r[io->nr].len;
		}
	}
#if SINGLE_THREADED
	iothread_doit(io);
#else
	io->state = WORK;
	sem_V(&io->start);
#endif
}

void iothread_read(iothread_t *io, int fd, void *buf, unsigned long off, unsigned long len, int *clen)
{
	iothread_start(io, fd, buf, 1, (int)off, (int)len, clen, NULL);
}

void iothread_write(iothread_t *io, int fd, void *buf, unsigned long off, unsigned long len, int *clen)
{
	iothread_start(io, fd, buf, 0, (int)off, (int)len, clen, NULL);
}

/* if the data has been compressed and was read, then this function must go
   through and uncompress each of the blocks that completed */

int iothread_wait(iothread_t *io, int check)
{
	int i;

#if !SINGLE_THREADED
	assert(io->state == WORK || io->state == REPLY);
	sem_P(&io->reply);
	assert(io->state == REPLY);
#endif
	if (check && io->res != io->total_len)
	{
		dump_iothread(io);
		if (io->res > 0)
			fprintf(stderr, "Incorrect io operation %ld != %ld\n", io->res, io->total_len);
		else if (!io->res)
			fprintf(stderr, "Unexpected eof.\n");
		else
			fprintf(stderr, "error #%ld\n", -io->res);
		exit(1);
	}

	for (i = 0; i < io->nr; i++)
	{
		if (compression_bits && io->r[i].is_read && io->r[i].clen)
		{
			int x, y;
			char *buf = io->r[i].buf;

			for (x = 0, y = io->r[i].len - *io->r[i].clen; y < io->r[i].len; x += BYTES1)
			{
				int compressed_bytes, l;

				buf[x] = buf[y++];
				compressed_bytes = ((buf[x] & ((1 << compression_bits) - 1)) * k - compression_bits - select_bits) / 8;
				if (!compressed_bytes)
					compressed_bytes = 1;
				for (l = compressed_bytes; l < BYTES1; l++)
					buf[x + l] = buf[y++];
			}
		}
	}

	io->state = START;
	return io->res;
}

void thd_init(void)
{
	int i;

	for (i = 0; i < 3; i++)
		iothread_create(&iot[i]);
}

int log2anf(unsigned long x)
{
	unsigned long y = 1;
	int i = 0;

	if (x & 0x80000000)
		return 32;
	while (y < x)
		y <<= 1, i++;
	return i;
}

void init_tables(FILE *startf, FILE *probef, int do_sim)
{
	unsigned long have, used;

	logn = log2anf(n);
	if (do_compression)
		compression_bits = log2anf(logn + bits);
	else
		compression_bits = 0;

	if (startf || probef)
		select_bits = 1;
	else
		select_bits = 0;

	//printf("logn = %d, bits = %d => total bits %d\n", logn, bits, logn+bits);

	assert(compression_bits <= 8);
	have = mem * 1024 * 1024;

#if EXPERIMENTAL_DIST
	if (do_sim)
	{
		dist = malloc(sizeof(*dist));
		have -= sizeof(*dist);
	}
#endif
	n_buckets[1] = 1;

	if (have <= BYTES(n))
	{
		for (;;)
		{
			int open_max = sysconf(_SC_OPEN_MAX);

			if (n_buckets[1] > 1)
				used = 2 * BYTES_IN(n_buckets[1]);
			else
				used = BYTES_IN(n_buckets[1]);
			n_buckets[0] = (BYTES(n) + have - 1 - used) / (have - used);
			if (n_buckets[0] > 1)
				n_buckets[0] *= 2;

			if (have >= used && n_buckets[0] * n_buckets[1] < open_max - 6)
				break;

			if (n_buckets[1] >= open_max - 6)
			{
				fprintf(stderr, "*** oops, not enough files are available for opening (need %d MB tried %d,%d).\n", BYTES(n) / (1024 * 1024), n_buckets[0], n_buckets[1]);
				exit(1);
			}

			n_buckets[1]++;
		}

		buffers[2] = malloc(BYTES_IN(n_buckets[1]));
		buffers[3] = malloc(BYTES_IN(n_buckets[1]));
		used = 2 * BYTES_IN(n_buckets[1]);
	}
	else
	{
		buffers[2] = malloc(BYTES(n));
		buffers[3] = NULL;
		used = BYTES(n);
	}

	have -= used;

	n_buckets[0] = (BYTES(n) + have - 1) / have; /* *2 to give me double buff */
	if (n_buckets[0] > 1)
		n_buckets[0] *= 2;

	if (verbose)
		fprintf(stderr, "Picked %d past = %d bytes, %d cur = %d bytes\n", n_buckets[0], BYTES_IN(n_buckets[0]), n_buckets[1], BYTES_IN(n_buckets[1]));

	/* make sure the buffers aren't too big or else we'll stall a lot on the
       last one of each iteration */
	if (n_buckets[0] > 1 && n_buckets[0] < 16)
		n_buckets[0] = 16;

	buffers[0] = malloc(BYTES_IN(n_buckets[0]));
	buffers[1] = NULL;
	have -= BYTES_IN(n_buckets[0]);

	if (n_buckets[0] > 1)
	{
		buffers[1] = malloc(BYTES_IN(n_buckets[0]));
		have -= BYTES_IN(n_buckets[0]);
		if (pastfname == NULL)
			pastfname = tempnam("/tmp/", "thd");
		pastfd = open(pastfname, O_RDWR | O_CREAT | O_TRUNC, 00700);
		unlink(pastfname);
		if (pastfd == -1)
		{
			perror(pastfname);
			exit(1);
		}
		/*
	lseek(pastfd, BYTES_IN(n_buckets[0])*n_buckets[0], SEEK_SET);
	write(pastfd, &zero, 1);
	close(pastfd);
	pastfd = open(pastfname, O_RDWR);
	assert(pastfd >= 0);
*/
	}
	if (n_buckets[1] > 1)
	{
		if (curfname == NULL)
			curfname = tempnam("/tmp/", "thd");
		curfd = open(curfname, O_RDWR | O_CREAT | O_TRUNC, 00700);
		unlink(curfname);
		if (curfd == -1)
		{
			perror(curfname);
			exit(1);
		}
	}

	past_clens = calloc(sizeof(*past_clens), n_buckets[0]);
	cur_clens = calloc(sizeof(*cur_clens), n_buckets[1]);

	edges_bytes_avail = have;
	BYTES1 = BYTES(1);

	if (compression_bits + select_bits > 8)
	{
		fprintf(stderr, "\n\n  **** WARNING ****\n\nThis code is currently untested for huge numbers of nodes.\nI suspect that the results may be incorrect.\n\n");
	}
}

static void
setup_predefined_edgefile(FILE *f)
{
	edgef = malloc(sizeof(*edgef));
	if (f)
		edgef[0] = f;
	else if ((edgef[0] = fopen(edgefname, "r+")) == NULL)
	{
		perror(edgefname);
		exit(1);
	}
	nedgef = 1;
	pastind = malloc(sizeof(*pastind));
	pastind[0] = 0;
	curind = pastind;
}

void read_graph(FILE *f, FILE *startf, FILE *probef, int do_dist)
{
	int i, j;
	char line[16 * 1024];
	edge_t edge;

	if (n <= 0)
	{
		if (fgets(line, sizeof(line), f) == NULL || sscanf(line, "%d", &n) != 1)
		{
			fprintf(stderr, "Error: first line of graph file must be # of nodes.\n");
			exit(1);
		}
	}

	if (n <= 0)
	{
		fprintf(stderr, "invalid number of nodes specified: %d\n", n);
		exit(1);
	}

	init_tables(startf, probef, do_dist);

	if (n_buckets[0] == 1)
	{
		/* there is some hope of fitting the edge file into memory */
		if (edgefname)
		{
			FILE *f = fopen(edgefname, "r");
			struct stat sbuf;

			if (!f)
			{
				perror(edgefname);
				exit(1);
			}

			if (stat(edgefname, &sbuf) == -1 || sbuf.st_size <= 0)
			{
				perror("stat");
				exit(1);
			}

			if (sbuf.st_size <= edges_bytes_avail)
			{
				edges = malloc(sbuf.st_size);
				total_edges = sbuf.st_size / sizeof(edges[0]);

				if (verbose)
					fprintf(stderr, "Loading disk based edge file into RAM.\n");
				if (fread(edges, 1, sbuf.st_size, f) != sbuf.st_size)
				{
					perror("fread");
					exit(1);
				}
				fclose(f);
				edges_bytes_avail -= sbuf.st_size;
			}
			else
			{
				if (verbose)
					fprintf(stderr, "Using existing disk based edge file.\n");
				setup_predefined_edgefile(f);
			}
		}
		else
		{
			/* start reading until we run out of memory at which point
	       create a temporary file and dump the edges that we've
	       already read */
			edges = malloc(edges_bytes_avail);
			total_edges = 0;

			if (verbose)
				fprintf(stderr, "Attempting to load edges into RAM, can take %d.\n", edges_bytes_avail / sizeof(*edges));

			while (fgets(line, sizeof(line), f) != NULL)
			{
				int u, v;

				if (sscanf(line, "%d %d", &u, &v) != 2)
				{
					fprintf(stderr, "Invalid input: %s", line);
					exit(1);
				}

				if (total_edges == edges_bytes_avail / sizeof(*edges))
				{
					char *fname = tempnam("/tmp", "thd");
					FILE *f;

					if (verbose)
						fprintf(stderr, "Could not load edges into RAM, dumping to disk and continuing.\n");
					if ((f = fopen(fname, "w+")) == NULL)
					{
						perror(fname);
						exit(1);
					}
					unlink(fname);
					setup_predefined_edgefile(f);
					if (fwrite(edges, sizeof(*edges), total_edges, f) != total_edges)
					{
						fprintf(stderr, "Write to edge file failed (out of disk space?).\n");
						exit(1);
					}
					free(edges);
					edges = NULL;
				}
				if (total_edges >= edges_bytes_avail / sizeof(*edges))
				{
					edge_t e;
					e.u = u;
					e.v = v;
					if (fwrite(&e, sizeof(e), 1, edgef[0]) != 1)
					{
						fprintf(stderr, "Write to edge file failed (out of disk space?).\n");
						exit(1);
					}
				}
				else
				{
					edges[total_edges].u = u;
					edges[total_edges].v = v;
				}
				total_edges++;
			}
			if (total_edges <= edges_bytes_avail / sizeof(*edges))
			{
				edges_bytes_avail -= total_edges * sizeof(*edges);
			}
		}
	}
	else if (n_buckets[0] * n_buckets[1] < sysconf(_SC_OPEN_MAX) - 6)
	{
		/* must be on disk when doing the splits */
		/* fortunately, we can do it in one pass */
		if (verbose)
			fprintf(stderr, "[%lu] Splitting edges into %dx%d\n", T, n_buckets[0], n_buckets[1]);

		edgef = malloc(sizeof(*edgef) * n_buckets[0] * n_buckets[1]);
		pastind = malloc(sizeof(*pastind) * n_buckets[0] * n_buckets[1]);
		curind = malloc(sizeof(*curind) * n_buckets[0] * n_buckets[1]);

		for (i = 0; i < n_buckets[1]; i++)
		{
			for (j = 0; j < n_buckets[0]; j++)
			{
				char *fname = tempnam("/tmp", "thd");

				edgef[i * n_buckets[0] + j] = fopen(fname, "w+");
				unlink(fname);
				if (edgef[i * n_buckets[0] + j] == NULL)
				{
					perror(fname);
					exit(1);
				}
			}
		}

		if (edgefname)
		{
			if ((f = fopen(edgefname, "r")) == NULL)
			{
				perror(edgefname);
				exit(1);
			}
		}

		for (;;)
		{
			if (edgefname)
			{
				if (fread(&edge, 1, sizeof(edge), f) != sizeof(edge))
					break;
			}
			else
			{
				if (fgets(line, sizeof(line), f) == NULL)
					break;
				if (sscanf(line, "%d %d", &edge.u, &edge.v) != 2)
				{
					fprintf(stderr, "Couldn't read edge on line: %s", line);
					exit(1);
				}
			}
			if (edge.v < 0 || edge.v >= n || edge.u < 0 || edge.u >= n)
			{
				fprintf(stderr, "invalid edge: %s", line);
				fprintf(stderr, "   [ one point is not in the range 0..%d]\n", n - 1);
				exit(1);
			}

			if (fwrite(&edge, 1, sizeof(edge), edgef[BUCKET(edge.u, n_buckets[1], n) * n_buckets[0] + BUCKET(edge.v, n_buckets[0], n)]) != sizeof(edge))
			{
				perror("fwrite");
				exit(1);
			}
		}

		nedgef = 0;
		for (i = 0; i < n_buckets[1]; i++)
		{
			for (j = 0; j < n_buckets[0]; j++)
			{
				if (ftell(edgef[i * n_buckets[0] + j]) > 0)
				{
					edgef[nedgef] = edgef[i * n_buckets[0] + j];
					pastind[nedgef] = N_IN(n_buckets[0]) * j;
					curind[nedgef] = N_IN(n_buckets[1]) * i;
					nedgef++;
				}
			}
		}

		if (edgefname)
			fclose(f);

		if (verbose)
			fprintf(stderr, "[%lu] Ended up with %d files.\n", T, nedgef);
	}
	else
	{
		/* must split it into too many buckets to do it in one pass */
		fprintf(stderr, "About to implement this code ... %d %d \n", n_buckets[0], n_buckets[1]);
		exit(1);
	}
}

void init_random_bits(FILE *start, FILE *probe)
{
	unsigned char *last, *cur, *lastcur, *curcur;
	int i, j, ki, last_start = -1, last_probe = -1;
	int bit, newbit;
	int index1, index2;
	int anylast, anylastcur;
	unsigned long probed_nodes = 0, starting_nodes = 0;

	last = buffers[1];
	cur = buffers[0];

	lastcur = buffers[3];
	curcur = buffers[2];

	if (verbose)
		fprintf(stderr, "[%lu] Initializing random bits + doing first copy.\n", T);

	if (n_buckets[0] > 1)
		lseek(pastfd, 0, SEEK_SET);
	if (n_buckets[1] > 1)
		lseek(curfd, 0, SEEK_SET);

	for (ki = 0; ki < nk; ki++)
	{
		anylast = 0;
		anylastcur = 0;

		index1 = 0;
		index2 = 0;

		memset(cur, 0, BYTES_IN(n_buckets[0]));
		memset(curcur, 0, BYTES_IN(n_buckets[1]));

		for (i = 0; i < n; i++)
		{
			bit = compression_bits + select_bits;
			if (probe)
			{
				while (last_probe < i)
				{
					int tmpi;

					if (fscanf(probe, "%d", &tmpi) != 1)
					{
						if (feof(probe))
							last_probe = n + 1;
						else
						{
							fprintf(stderr, "Error reading from probe set at node %d\n", i);
							exit(1);
						}
					}
					else
					{
						if (tmpi < last_probe)
						{
							fprintf(stderr, "Error: probe set is not sorted!\n");
							exit(1);
						}
						last_probe = tmpi;
					}
				}
				if (!start)
					last_start = last_probe;
			}
			if (start)
			{
				while (last_start < i)
				{
					int tmpi;

					if (fscanf(start, "%d", &tmpi) != 1)
					{
						if (feof(start))
							last_start = n + 1;
						else
						{
							fprintf(stderr, "Error reading from start set at node %d\n", i);
							exit(1);
						}
					}
					else
					{
						if (tmpi < last_start)
						{
							fprintf(stderr, "Error: start set is not sorted!\n");
							exit(1);
						}
						last_start = tmpi;
					}
				}
				if (!probe)
					last_probe = last_start;
			}

			for (j = 0; j < k; j++)
			{
				unsigned long x;
				int actual_bit;

				if (j == 0 && select_bits)
				{
					if (last_start == i)
					{
						cur[index1 + compression_bits / 8] |= (1 << (compression_bits % 8));
						curcur[index2 + compression_bits / 8] |= (1 << (compression_bits % 8));
						starting_nodes++;
					}
					if (last_probe != i)
					{
						break;
					}
				}

				if (j == 0)
					probed_nodes++;

				x = random();
#if STOCHASTIC
				j = random() % k;
#endif
				for (newbit = 0; (x & (1 << newbit)) == 0 && newbit < logn + bits; newbit++)
				{
				}
				if (newbit >= logn + bits)
					newbit = logn + bits - 1;

				if (compression_bits)
				{
					actual_bit = newbit * k + j + select_bits + compression_bits;
				}
				else
				{
					actual_bit = newbit + select_bits + j * (logn + bits);
				}

				cur[index1 + actual_bit / 8] |= (1 << (actual_bit % 8));
				curcur[index2 + actual_bit / 8] |= (1 << (actual_bit % 8));
#if STOCHASTIC
				break;
#endif
			}

			index1 += BYTES1;
			index2 += BYTES1;

			if (n_buckets[0] > 1 && (i + 1) % N_IN(n_buckets[0]) == 0)
			{
				unsigned char *tmp = last;

				if (anylast)
					iothread_wait(&iot[0], 1);
				anylast = 1;
				iothread_write(&iot[0], pastfd, cur, FILE_OFFSET(i, n_buckets[0], n), BYTES_IN(n_buckets[0]), &past_clens[BUCKET(i, n_buckets[0], n)]);
				last = cur;
				cur = tmp;
				index1 = 0;
				memset(cur, 0, BYTES_IN(n_buckets[0]));
			}
			if (n_buckets[1] > 1 && (i + 1) % N_IN(n_buckets[1]) == 0)
			{
				unsigned char *tmp = lastcur;
				if (anylastcur)
					iothread_wait(&iot[1], 1);
				anylastcur = 1;
				iothread_write(&iot[1], curfd, curcur, FILE_OFFSET(i, n_buckets[1], n), BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(i, n_buckets[1], n)]);
				lastcur = curcur;
				curcur = tmp;
				index2 = 0;
				memset(curcur, 0, BYTES_IN(n_buckets[1]));
			}
		}
		if (n_buckets[0] > 1 && index1)
		{
			if (anylast)
				iothread_wait(&iot[0], 1);
			iothread_write(&iot[0], pastfd, cur, FILE_OFFSET(i, n_buckets[0], n), BYTES_IN(n_buckets[0]), &past_clens[BUCKET(i, n_buckets[0], n)]);
			anylast = 1;
		}
		if (n_buckets[1] > 1 && index2)
		{
			if (anylastcur)
				iothread_wait(&iot[1], 1);
			iothread_write(&iot[1], curfd, curcur, FILE_OFFSET(i, n_buckets[1], n), BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(i, n_buckets[1], n)]);
			anylastcur = 1;
		}
		if (anylast)
			iothread_wait(&iot[0], 1);
		if (anylastcur)
			iothread_wait(&iot[1], 1);
	}
	if ((starting_nodes || probed_nodes) && verbose)
		fprintf(stderr, "From %ld probing %ld of %d nodes.\n", starting_nodes, probed_nodes, n);
#if EXPERIMENTAL_DIST
	if (dist)
	{
		int nodes = starting_nodes;
		if (!nodes)
			nodes = n;
		anf_dist_init(dist, nodes, edges_bytes_avail);
	}
#endif
}

static void
fix_compressed_bits(unsigned char *cur, int n_buckets)
{
	if (compression_bits)
	{
		int j;
		int max = BYTES_IN(n_buckets);
		int bytes = BYTES1;

		for (j = 0; j < max; j += bytes, cur += bytes)
		{
			int compressed_ones = cur[0] & ((1 << compression_bits) - 1);
			int first_bit = compression_bits + compressed_ones * k + select_bits;
			int i;

			for (i = first_bit;; i += k, compressed_ones++)
			{
				int j;

				for (j = 0; j < k && (cur[(i + j) / 8] & (1 << ((i + j) % 8))); j++)
				{
				}
				if (j < k)
					break;
			}
			if (i > first_bit)
			{
				cur[0] = (cur[0] & ~((1 << compression_bits) - 1)) | compressed_ones;
			}
		}
	}
}

void do_iteration(void)
{
	edge_t edge;
	int i, ki = 0, ef;
	int nedges = 0;
	int pre;
	unsigned char *cur, *last, *tmp;
	unsigned char *curcur, *readwritecur;
	static int it = 0;
	int compressed_ones;
	int first_bit;

	if (!edges && !nedgef)
		return;

	it++;

	if (verbose)
		fprintf(stderr, "[%lu] Starting iteration %d.%d\n", T, it, ki);
	if (n_buckets[0] == 1)
	{
		int bytes = BYTES1;

		assert(n_buckets[1] == 1);
		if (!edges)
			rewind(edgef[0]);
		for (;;)
		{
			if (edges)
			{
				if (nedges >= total_edges)
					break;
				cur = &buffers[2][bytes * edges[nedges].u];
				last = &buffers[0][bytes * edges[nedges].v];
			}
			else
			{
				if (fread(&edge, 1, sizeof(edge), edgef[0]) <= 0)
					break;
				cur = &buffers[2][bytes * edge.u];
				last = &buffers[0][bytes * edge.v];
			}

			/* n.b. we can assume that cur & last have same alignment
	           because malloc returns things aligned and they are in sync */

			nedges++;
#if 1
			if (compression_bits)
			{
				int x = last[0] & ((1 << compression_bits) - 1);
				compressed_ones = cur[0] & ((1 << compression_bits) - 1);
				if (x > compressed_ones)
				{
					compressed_ones = x;
					cur[0] = (cur[0] & ~((1 << compression_bits) - 1)) | compressed_ones;
				}
				first_bit = select_bits + compression_bits + compressed_ones * k;
			}
			else
			{
				first_bit = select_bits;
			}

			i = first_bit / 8;
			if (i == 0 && first_bit)
			{
				cur[0] |= last[0] & ~((1 << (compression_bits + select_bits)) - 1);
				i = 1;
			}

			pre = 8 - (((unsigned long)&cur[i]) & 0x7);
			if (pre == 8)
				pre = 0;
			pre += i;

			for (; i < bytes && i < pre; i++)
				cur[i] |= last[i];
			while (i + 4 <= bytes)
			{
				(*(unsigned long *)&cur[i]) |= (*(unsigned long *)&last[i]);
				i += 4;
			}
			while (i < bytes)
			{
				cur[i] |= last[i];
				i++;
			}
#else
			for (i = 0; i < bytes; i += 4)
				(*(unsigned long *)&cur[i]) |= (*(unsigned long *)&last[i]);
#endif
		}

		fix_compressed_bits(buffers[2], n_buckets[1]);
	}
	else
	{
		cur = buffers[0];
		last = buffers[1];

		curcur = buffers[3];
		readwritecur = buffers[2];

		iothread_read(&iot[0], pastfd, last, BYTES_IN(n_buckets[0]) * BUCKET(pastind[0], n_buckets[1], n), BYTES_IN(n_buckets[0]), &past_clens[BUCKET(pastind[0], n_buckets[0], n)]);
		if (n_buckets[1] > 1)
		{
			iothread_read(&iot[1], curfd, readwritecur, BYTES_IN(n_buckets[0]) * BUCKET(curind[0], n_buckets[1], n), BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[0], n_buckets[1], n)]);
		}

		for (ef = 0; ef < nedgef; ef++)
		{
			rewind(edgef[ef]);

			if (ef == 0 || pastind[ef - 1] != pastind[ef])
			{
				tmp = cur;
				cur = last;
				last = tmp;
				iothread_wait(&iot[0], 1);
				for (i = ef + 1; i < nedgef; i++)
				{
					if (pastind[i - 1] != pastind[i])
					{
						iothread_read(&iot[0], pastfd, last, BYTES1 * pastind[i], BYTES_IN(n_buckets[0]), &past_clens[BUCKET(pastind[i], n_buckets[0], n)]);
						break;
					}
				}
			}

			if (ef == 0 || curind[ef - 1] != curind[ef])
			{
				tmp = curcur;
				curcur = readwritecur;
				readwritecur = tmp;

				if (n_buckets[1] > 1)
					iothread_wait(&iot[1], 1);

				for (i = ef + 1; i < nedgef; i++)
				{
					if (curind[i - 1] != curind[i])
					{
						if (ef)
						{
							fix_compressed_bits(readwritecur, n_buckets[1]);
							iothread_start(&iot[1], curfd, readwritecur, 0, BYTES1 * curind[ef - 1], BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[ef - 1], n_buckets[1], n)], readwritecur, 1, BYTES1 * curind[i], BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[i], n_buckets[1], n)], NULL);
						}
						else
						{
							iothread_read(&iot[1], curfd, readwritecur, BYTES1 * curind[i], BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[i], n_buckets[1], n)]);
						}
						break;
					}
				}
				if (i >= nedgef && n_buckets[1] > 1)
				{
					fix_compressed_bits(readwritecur, n_buckets[1]);
					iothread_write(&iot[1], curfd, readwritecur, BYTES1 * curind[ef - 1], BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[ef - 1], n_buckets[1], n)]);
				}
			}

			while (fread(&edge, 1, sizeof(edge), edgef[ef]) > 0)
			{
				unsigned char *curp = &curcur[BYTES1 * (edge.u - curind[ef])];
				unsigned char *last = &cur[BYTES1 * (edge.v - pastind[ef])];
				int bytes = BYTES1;

				assert(edge.u >= curind[ef]);
				assert(edge.u < curind[ef] + N_IN(n_buckets[1]));
				assert(edge.v >= pastind[ef]);
				assert(edge.v < pastind[ef] + N_IN(n_buckets[0]));

				nedges++;
#if 1
				if (compression_bits)
				{
					int x = last[0] & ((1 << compression_bits) - 1);
					compressed_ones = curp[0] & ((1 << compression_bits) - 1);
					if (x > compressed_ones)
					{
						compressed_ones = x;
						curp[0] = (curp[0] & ~((1 << compression_bits) - 1)) | compressed_ones;
					}
					first_bit = select_bits + compression_bits + compressed_ones * k;
				}
				else
				{
					first_bit = select_bits;
				}

				i = first_bit / 8;
				if (i == 0 && first_bit)
				{
					curp[0] |= last[0] & ~((1 << (select_bits + compression_bits)) - 1);
					i = 1;
				}

				pre = 8 - (((unsigned long)&curp[i]) & 0x7);
				if (pre == 8)
					pre = 0;
				pre += i;

				for (; i < bytes && i < pre; i++)
					curp[i] |= last[i];
				while (i + 4 <= bytes)
				{
					(*(unsigned long *)&curp[i]) |= (*(unsigned long *)&last[i]);
					i += 4;
				}
				while (i < bytes)
				{
					curp[i] |= last[i];
					i++;
				}
#else
				for (i = 0; i < bytes; i += 4)
					(*(unsigned long *)&curp[i]) |= (*(unsigned long *)&last[i]);
#endif
			}
		}
		fix_compressed_bits(curcur, n_buckets[1]);
		if (n_buckets[1] > 1)
		{
			iothread_wait(&iot[1], 1);
			iothread_write(&iot[1], curfd, curcur, BYTES1 * curind[nedgef - 1], BYTES_IN(n_buckets[1]), &cur_clens[BUCKET(curind[nedgef - 1], n_buckets[1], n)]);
			iothread_wait(&iot[1], 1);
		}
	}
}

double
make_estimate(double *full_est)
{
	int i, j, b;
	double est = 0;
	unsigned char *cur0, *cur = buffers[2], *next = buffers[3];
	unsigned char *past = buffers[0];
	unsigned char *tmp;
	int bit;
	double full = 0;
	double tmppow;
	int dist_point = 0;

	if (verbose)
		fprintf(stderr, "[%lu] Computing estimate.\n", T);

	if (pastfd >= 0)
		lseek(pastfd, 0, SEEK_SET);
	if (curfd >= 0)
		lseek(curfd, 0, SEEK_SET);

	if (n_buckets[0] == 1)
		memcpy(buffers[0], cur, BYTES1 * n);
	else if (n_buckets[1] == 1)
	{
		if (!compression_bits)
		{
			iothread_write(&iot[1], pastfd, cur, 0, BYTES_IN(n_buckets[1]), NULL);
		}
	}
	else
		iothread_read(&iot[0], curfd, next, 0, BYTES_IN(n_buckets[1]), &cur_clens[0]);

	for (b = 0; b < n_buckets[1]; b++)
	{
		int ind0 = N_IN(n_buckets[1]) * b;
		int indn = ind0 + N_IN(n_buckets[1]);

		if (indn > n)
			indn = n;
		if (b && !compression_bits)
			iothread_wait(&iot[1], 1);
		if (n_buckets[1] > 1)
		{
			iothread_wait(&iot[0], 1);
			tmp = cur;
			cur = next;
			next = tmp;
			if (!compression_bits)
				iothread_write(&iot[1], pastfd, cur, -1, BYTES_IN(n_buckets[1]), NULL);
			if (b < n_buckets[1] - 1)
				iothread_read(&iot[0], curfd, next, BYTES1 * N_IN(n_buckets[1]) * (b + 1), BYTES_IN(n_buckets[1]), &cur_clens[b + 1]);
		}

		cur0 = cur;

		for (i = ind0; i < indn; i++, cur0 += BYTES1)
		{
			double b = 0;
			int compressed_bits = 0;

			if (compression_bits)
			{
				compressed_bits = (cur0[0] & ((1 << compression_bits) - 1));
				if (n_buckets[0] > 1)
				{
					int pastind = i % N_IN(n_buckets[0]);
					wordcopy(&past[BYTES1 * pastind], cur0, BYTES1);
					if (pastind == N_IN(n_buckets[0]) - 1 || i == n - 1)
					{
						if (pastind != i)
							iothread_wait(&iot[1], 1);
						iothread_write(&iot[1], pastfd, past, FILE_OFFSET(i, n_buckets[0], n), BYTES_IN(n_buckets[0]), &past_clens[BUCKET(i, n_buckets[0], n)]);
						if (past == buffers[0])
							past = buffers[1];
						else
							past = buffers[0];
					}
				}
			}

			tmppow = 0;
			for (j = 0; j < k; j++)
			{
				bit = compressed_bits;
				while (bit < logn + bits)
				{
					int want_bit;

					if (compression_bits)
						want_bit = compression_bits + bit * k + j;
					else
						want_bit = j * (logn + bits) + bit;
					want_bit += select_bits;

					if ((cur0[want_bit / 8] & (1 << (want_bit % 8))) == 0)
					{
						break;
					}
					bit++;
				}
				if (bit >= logn + bits)
					bit = (logn + bits);
				b += bit;
			}

#if STOCHASTIC
			tmppow = (int)((k / (.77351 * BIAS)) * pow(2, b / k));
#else
			tmppow = pow(2.0, b / k);
#endif
			if (!select_bits || (cur0[compression_bits / 8] & (1 << compression_bits)) != 0)
			{
				est += tmppow;
#if EXPERIMENTAL_DIST
				if (dist)
					anf_dist_point(dist, dist_point++, tmppow);
#endif
			}
			full += tmppow;
		}
	}
#if EXPERIMENTAL_DIST
	if (dist)
		anf_dist_iteration_done(dist);
#endif
	if (n_buckets[0] > 1)
		iothread_wait(&iot[1], 1);
#if STOCHASTIC
	if (full_est)
		*full_est = full;
#else
	if (full_est)
		*full_est = full / .77351;
#endif
#if STOCHASTIC
	return est;
#else
	return est / .77351;
#endif
}

double
individual_estimate(int node)
{
	int j;
	double b = 0;
	unsigned char *cur0, *cur = buffers[2];
	int bit;

	if (n_buckets[0] != 1)
	{
		fprintf(stderr, "Sorry, this only works for in-core processing.\n");
		exit(1);
	}

	if (compression_bits)
	{
		fprintf(stderr, "Sorry, this only works for uncompressed processing.\n");
		exit(1);
	}

	cur0 = &cur[BYTES(node)];
	if (select_bits && (cur0[compression_bits / 8] & (1 << (compression_bits % 8))) == 0)
	{
		return 0;
	}

	for (j = 0; j < k; j++)
	{
		bit = j * (logn + bits) + select_bits;
		while (bit < (j + 1) * (logn + bits) + select_bits)
		{
			if ((cur0[bit / 8] & (1 << (bit % 8))) == 0)
			{
				b += (bit - select_bits) % (logn + bits);
				break;
			}
			bit++;
		}
		if (bit == (j + 1) * (logn + bits) + select_bits)
			b += (logn + bits);
	}

	return pow(2.0, b / k) / (.7731 * BIAS);
}

double
jaccard(int u, int v)
{
	int j;
	double b = 0;
	unsigned char *u0, *v0, *cur = buffers[2];
	int bit;
	double size_u = individual_estimate(u);
	double size_v = individual_estimate(v);
	double size_union;
	double size_intersect;

	if (n_buckets[0] != 1)
	{
		fprintf(stderr, "Sorry, this only works for in-core processing.\n");
		exit(1);
	}

	if (compression_bits)
	{
		fprintf(stderr, "Sorry, this only works for uncompressed processing.\n");
		exit(1);
	}

	u0 = &cur[BYTES(u)];
	v0 = &cur[BYTES(v)];

	for (j = 0; j < k; j++)
	{
		bit = j * (logn + bits) + select_bits;
		while (bit < (j + 1) * (logn + bits) + select_bits)
		{
			if ((u0[bit / 8] & (1 << (bit % 8))) == 0 && (v0[bit / 8] & (1 << (bit % 8))) == 0)
			{
				b += (bit - select_bits) % (logn + bits);
				break;
			}
			bit++;
		}
		if (bit == (j + 1) * (logn + bits) + select_bits)
			b += (logn + bits);
	}

	size_union = pow(2.0, b / k) / (.7731 * BIAS);
	size_intersect = size_u + size_v - size_union;

	return size_intersect / size_union;
}

#if 0
void
check()
{
     int count[logn+bits];
     int i, j, kk;
     unsigned char *cur, *last, *tmp;

     if (verbose) fprintf(stderr, "[%lu] Starting counting.\n", T);

     for (i = 0; i < logn+bits; i++) count[i] = 0;
     cur = buffers[1];
     last = buffers[0];
     if (pastfd >= 0) lseek(pastfd, 0, SEEK_SET);
     if (n_buckets[0] > 1) reader_start(pastfd, last, -1, BYTES_IN(n_buckets[0]));

     for (i = 0; i < n_buckets[0]; i++) {
	if (n_buckets[0] > 1) reader_wait(1);
	if (i < n_buckets[0]-1) {
	    tmp  = last;
	    last = cur;
	    cur  = tmp;
	    reader_start(pastfd, last, -1, BYTES_IN(n_buckets[0]));
	} else {
	    cur = last;
	}
	for (j = 0; j < N_IN(n_buckets[0]); j++) {
	    for (kk = 0; kk < k*(logn+bits); kk++) {
		if (cur[j*BYTES1+kk/8] & (1 << (kk%8))) {
		    count[kk%(logn+bits)]++;
		}
	    }
	}
    }

    if (verbose) fprintf(stderr, "[%lu] Done counting.\n", T);

    for (i = 0; i < logn+bits; i++) {
	if (count[i]) printf("%d %d\n", i, count[i]);
    }

     if (verbose) fprintf(stderr, "[%lu] Starting counting.\n", T);

     for (i = 0; i < logn+bits; i++) count[i] = 0;
     cur = buffers[3];
     last = buffers[2];
     if (curfd >= 0) lseek(curfd, 0, SEEK_SET);
     if (n_buckets[1] > 1) reader_start(curfd, last, -1, BYTES_IN(n_buckets[1]));

     for (i = 0; i < n_buckets[1]; i++) {
	if (n_buckets[1] > 1) reader_wait(1);
	if (i < n_buckets[1]-1) {
	    tmp  = last;
	    last = cur;
	    cur  = tmp;
	    reader_start(curfd, last, -1, BYTES_IN(n_buckets[1]));
	} else {
	    cur = last;
	}
	for (j = 0; j < N_IN(n_buckets[1]); j++) {
	    for (kk = 0; kk < k*(logn+bits); kk++) {
		if (cur[j*BYTES1+kk/8] & (1 << (kk%8))) {
		    count[kk%(logn+bits)]++;
		}
	    }
	}
    }

    if (verbose) fprintf(stderr, "[%lu] Done counting.\n", T);

    for (i = 0; i < logn+bits; i++) {
	if (count[i]) printf("%d %d\n", i, count[i]);
    }
}

#define N (1024 * 1024)

void
copyfiles(FILE *fr, FILE *fw)	/* copy fr to fw */
{
    int nb;
    char *tmp;
    unsigned char junk[N], junk1[N];
    char *a, *b;

    rewind(fr);
    rewind(fw);

#if 1
    a = junk;
    b = junk1;
    reader_start(fr, a, 0, N);
    while ((nb = reader_wait(0)) > 0) {
	writer_start(fw, a, -1, nb);
	reader_start(fr, b, -1, N);
	writer_wait(1);
	tmp = b;
	a = b;
	b = tmp;
    }
#else
    while ((nb = fread(junk, 1, N, fr)) > 0) {
	fwrite(junk, 1, nb, fw);
    }
#endif
}

#endif

double *estimates[10 * 1024];
#define EFF_DIAM_FRAC .9

int main(int argc, char **argv)
{
	int it, maxit = -1;
	double lastest = -1, est = 0, full_est = 0;
	unsigned long seed;
	int individual = 0;
	int similarities = 0;
	int distances = 0;
	int sim_u = 0, sim_v = 0;
	int hog_mem = 0;
	FILE *startf, *probef;
	char *startfname = NULL, *probefname = NULL;

	seed = T;

	while (argc > 1 && argv[1][0] == '-')
	{
		if (argc > 2 && strcmp(argv[1], "-mem") == 0 && sscanf(argv[2], "%d", &mem) == 1)
			argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-hog") == 0 && sscanf(argv[2], "%d", &hog_mem) == 1)
			argc--, argv++;
		else if (strcmp(argv[1], "-verbose") == 0)
			verbose = 1;
		else if (argc > 2 && strcmp(argv[1], "-k") == 0 && sscanf(argv[2], "%d", &k) == 1)
			nk = 1, argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-n") == 0 && sscanf(argv[2], "%d", &n) == 1)
			argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-past") == 0)
			pastfname = argv[2], argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-cur") == 0)
			curfname = argv[2], argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-edges") == 0)
			edgefname = argv[2], argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-probe") == 0)
			probefname = argv[2], argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-start") == 0)
			startfname = argv[2], argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-bits") == 0 && sscanf(argv[2], "%d", &bits) == 1)
			argc--, argv++;
		else if (argc > 2 && strcmp(argv[1], "-it") == 0 && sscanf(argv[2], "%d", &maxit) == 1)
			argc--, argv++;
#if EXPERIMENTAL_DIST
		else if (strcmp(argv[1], "-dist") == 0)
			distances = 1;
#endif
		else if (strcmp(argv[1], "-compress") == 0)
			do_compression = 1;
		else if (strcmp(argv[1], "-norandom") == 0)
			seed = 0;
		else if (strcmp(argv[1], "-individual") == 0)
			individual = 1;
		else if (argc > 3 && strcmp(argv[1], "-sim") == 0 && sscanf(argv[2], "%d", &sim_u) == 1 && sscanf(argv[3], "%d", &sim_v) == 1)
		{
			similarities = 1;
			argc -= 2;
			argv += 2;
		}
		else
		{
			fprintf(stderr, "usage: [-mem MB | -k k | -bits bits | -it it | -past fname | -n n | -edges edgesfile | -past pastfname | -cur curfname  | -individual | -sim u v | -verbose | -compress | -hog MB | -start sorted-node-list-file | -probe sorted-node-list]\n");
			exit(1);
		}
		argc--, argv++;
	}

	if (similarities + individual + distances > 1)
	{
		fprintf(stderr, "May only specify one of similarities, individual and distance computations.\n");
		exit(1);
	}

	if (hog_mem > 0)
	{
		char *buf = malloc(hog_mem * 1024 * 1024);
		int i;
		for (i = 0; i < hog_mem * 1024 * 1024; i++)
			buf[i] = i;
	}

	if (startfname)
	{
		if ((startf = fopen(startfname, "r")) == 0)
		{
			perror(startfname);
			exit(1);
		}
	}
	else
	{
		startf = NULL;
	}

	if (probefname)
	{
		if ((probef = fopen(probefname, "r")) == 0)
		{
			perror(probefname);
			exit(1);
		}
	}
	else
	{
		probef = NULL;
	}

	srandom(seed);

	assert(sizeof(unsigned long) == 8);

	read_graph(stdin, startf, probef, distances);
	thd_init();
	init_random_bits(startf, probef);

	for (it = 1; (maxit <= 0 || it <= maxit) && (maxit > 0 || lastest * 1.01 < full_est); it++)
	{
		if (compression_bits && n_buckets[0] > 1 && verbose)
		{
			int i, len = 0;
			int total = BYTES_IN(n_buckets[0]) * n_buckets[0];

			if (n_buckets[1] > 1)
				total += BYTES_IN(n_buckets[1]) * n_buckets[1];

			for (i = 0; i < n_buckets[0]; i++)
				len += past_clens[i];
			for (i = 0; i < n_buckets[1]; i++)
				len += cur_clens[i];
			fprintf(stderr, "current size %.2f%%\n", ((double)len) / total * 100);
		}
		do_iteration();
		lastest = full_est;
		if (similarities)
		{
			printf("%d %f\n", it, jaccard(sim_u, sim_v));
			est = make_estimate(&full_est);
		}
		else if (!individual)
		{
			est = make_estimate(&full_est);
			printf("%d %f\n", it, est);
			fflush(stdout);
		}
		else
		{
			int i;

			estimates[it] = malloc(sizeof(*estimates[it]) * n);
			for (i = 0; i < n; i++)
			{
				estimates[it][i] = individual_estimate(i);
			}
			est = make_estimate(&full_est);
		}
	}

#if EXPERIMENTAL_DIST
	if (distances)
		anf_dist_printf(dist, stdout);
#endif

	if (individual)
	{
		int i, j;

		for (i = 0; i < n; i++)
		{
			int effdiam;
			double x_square = 0, x_sum = 0, x_y = 0, y_sum = 0;

			if (estimates[1][i] == 0)
			{
				/* we didn't do this node */
				continue;
			}

			for (j = 1; estimates[j][i] < estimates[it - 1][i] * EFF_DIAM_FRAC; j++)
			{
			}
			effdiam = j;
			for (j = 1; j <= effdiam; j++)
			{
				double x = log(j);
				double y = log(estimates[j][i]);

				x_square += x * x;
				x_sum += x;
				x_y += x * y;
				y_sum += y;
			}
			if (effdiam == 1)
				printf("1 0 %d 0\n", i);
			else
				printf("%d %f %d %f\n", effdiam, (n * x_y - x_sum * y_sum) / (n * x_square - x_sum * x_sum), i, ((double)estimates[effdiam][i]) / effdiam);
		}
	}
	return 0;
}
