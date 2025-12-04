/*******************************************************************************************
 *
 *  GIXpack: Compress/decompress ktab files using zstd seekable format
 *
 *  Author:  Erik Garrison (based on Gene Myers' FASTGA code)
 *  Date  :  December 2024
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "zstd.h"
#include "zstd_seekable.h"
#include "gene_core.h"

static char *Usage = "[-d] [-l<int(3)>] [-f<int(262144)>] <gix_path>";

// Compress a single ktab file to .zst format
// Format: [12-byte header uncompressed] [zstd seekable compressed data]
static int compress_ktab(char *inpath, char *outpath, int level, unsigned frameSize)
{
  FILE *fin, *fout;
  ZSTD_seekable_CStream *cstream;
  size_t const buffInSize = ZSTD_CStreamInSize();
  size_t const buffOutSize = ZSTD_CStreamOutSize();
  void *buffIn, *buffOut;
  char header[12];
  size_t toRead, ret;

  fin = fopen(inpath, "rb");
  if (fin == NULL)
    { fprintf(stderr, "Cannot open %s for reading\n", inpath);
      return 1;
    }

  fout = fopen(outpath, "wb");
  if (fout == NULL)
    { fprintf(stderr, "Cannot open %s for writing\n", outpath);
      fclose(fin);
      return 1;
    }

  // Read and write header uncompressed (4 bytes kmer + 8 bytes nels)
  if (fread(header, 1, 12, fin) != 12)
    { fprintf(stderr, "Cannot read header from %s\n", inpath);
      fclose(fin);
      fclose(fout);
      return 1;
    }
  if (fwrite(header, 1, 12, fout) != 12)
    { fprintf(stderr, "Cannot write header to %s\n", outpath);
      fclose(fin);
      fclose(fout);
      return 1;
    }

  // Initialize seekable compression
  cstream = ZSTD_seekable_createCStream();
  if (cstream == NULL)
    { fprintf(stderr, "Cannot create compression stream\n");
      fclose(fin);
      fclose(fout);
      return 1;
    }

  ret = ZSTD_seekable_initCStream(cstream, level, 0, frameSize);
  if (ZSTD_isError(ret))
    { fprintf(stderr, "Cannot init compression stream: %s\n", ZSTD_getErrorName(ret));
      ZSTD_seekable_freeCStream(cstream);
      fclose(fin);
      fclose(fout);
      return 1;
    }

  buffIn = malloc(buffInSize);
  buffOut = malloc(buffOutSize);
  if (buffIn == NULL || buffOut == NULL)
    { fprintf(stderr, "Cannot allocate buffers\n");
      ZSTD_seekable_freeCStream(cstream);
      fclose(fin);
      fclose(fout);
      free(buffIn);
      free(buffOut);
      return 1;
    }

  // Compress the data
  while ((toRead = fread(buffIn, 1, buffInSize, fin)) > 0)
    { ZSTD_inBuffer input = { buffIn, toRead, 0 };
      while (input.pos < input.size)
        { ZSTD_outBuffer output = { buffOut, buffOutSize, 0 };
          ret = ZSTD_seekable_compressStream(cstream, &output, &input);
          if (ZSTD_isError(ret))
            { fprintf(stderr, "Compression error: %s\n", ZSTD_getErrorName(ret));
              goto cleanup;
            }
          fwrite(buffOut, 1, output.pos, fout);
        }
    }

  // End the stream and write seek table
  { ZSTD_outBuffer output = { buffOut, buffOutSize, 0 };
    ret = ZSTD_seekable_endStream(cstream, &output);
    while (ret > 0)
      { fwrite(buffOut, 1, output.pos, fout);
        output.pos = 0;
        ret = ZSTD_seekable_endStream(cstream, &output);
      }
    fwrite(buffOut, 1, output.pos, fout);
  }

  free(buffIn);
  free(buffOut);
  ZSTD_seekable_freeCStream(cstream);
  fclose(fin);
  fclose(fout);
  return 0;

cleanup:
  free(buffIn);
  free(buffOut);
  ZSTD_seekable_freeCStream(cstream);
  fclose(fin);
  fclose(fout);
  return 1;
}

// Decompress a .zst ktab file back to original
static int decompress_ktab(char *inpath, char *outpath)
{
  FILE *fin, *fout;
  ZSTD_seekable *seekable;
  char header[12];
  int kmer;
  int64_t nels;
  (void)kmer; (void)nels;  // Used only for reading header
  size_t offset;
  void *buffer;
  size_t const buffSize = 1024 * 1024;  // 1MB chunks

  fin = fopen(inpath, "rb");
  if (fin == NULL)
    { fprintf(stderr, "Cannot open %s for reading\n", inpath);
      return 1;
    }

  // Read header
  if (fread(header, 1, 12, fin) != 12)
    { fprintf(stderr, "Cannot read header from %s\n", inpath);
      fclose(fin);
      return 1;
    }
  memcpy(&kmer, header, 4);
  memcpy(&nels, header + 4, 8);

  // Get file size to know compressed data size
  fseek(fin, 0, SEEK_END);
  long fileSize = ftell(fin);
  size_t compressedSize = fileSize - 12;
  fseek(fin, 12, SEEK_SET);

  // Read compressed data into memory for seekable init
  void *compressedData = malloc(compressedSize);
  if (compressedData == NULL)
    { fprintf(stderr, "Cannot allocate memory for compressed data\n");
      fclose(fin);
      return 1;
    }
  if (fread(compressedData, 1, compressedSize, fin) != compressedSize)
    { fprintf(stderr, "Cannot read compressed data\n");
      free(compressedData);
      fclose(fin);
      return 1;
    }
  fclose(fin);

  // Initialize seekable decompressor
  seekable = ZSTD_seekable_create();
  if (seekable == NULL)
    { fprintf(stderr, "Cannot create seekable decompressor\n");
      free(compressedData);
      return 1;
    }

  size_t ret = ZSTD_seekable_initBuff(seekable, compressedData, compressedSize);
  if (ZSTD_isError(ret))
    { fprintf(stderr, "Cannot init seekable: %s\n", ZSTD_getErrorName(ret));
      ZSTD_seekable_free(seekable);
      free(compressedData);
      return 1;
    }

  // Calculate data size from header (we need to know pbyte to calculate)
  // For now, decompress frame by frame until exhausted
  fout = fopen(outpath, "wb");
  if (fout == NULL)
    { fprintf(stderr, "Cannot open %s for writing\n", outpath);
      ZSTD_seekable_free(seekable);
      free(compressedData);
      return 1;
    }

  // Write header
  fwrite(header, 1, 12, fout);

  // Decompress all data
  buffer = malloc(buffSize);
  if (buffer == NULL)
    { fprintf(stderr, "Cannot allocate buffer\n");
      fclose(fout);
      ZSTD_seekable_free(seekable);
      free(compressedData);
      return 1;
    }

  offset = 0;
  while (1)
    { ret = ZSTD_seekable_decompress(seekable, buffer, buffSize, offset);
      if (ZSTD_isError(ret))
        { fprintf(stderr, "Decompression error: %s\n", ZSTD_getErrorName(ret));
          break;
        }
      if (ret == 0)
        break;
      fwrite(buffer, 1, ret, fout);
      offset += ret;
    }

  free(buffer);
  fclose(fout);
  ZSTD_seekable_free(seekable);
  free(compressedData);
  return 0;
}

int main(int argc, char *argv[])
{
  int decompress = 0;
  int level = 3;
  unsigned frameSize = 256 * 1024;  // 256KB frames
  char *gixpath;
  char *dir, *root;
  int p, nthreads;
  int f;

  // Parse arguments
  { int i, j;
    char *eptr;
    int flags[128];

    ARG_INIT("GIXpack")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { case 'd':
            decompress = 1;
            break;
          case 'l':
            ARG_NON_NEGATIVE(level, "l")
            if (level > 19)
              level = 19;
            break;
          case 'f':
            ARG_NON_NEGATIVE(frameSize, "f")
            if (frameSize < 1024)
              frameSize = 1024;
            break;
          default:
            fprintf(stderr, "%s: Unknown option %s\n", Prog_Name, argv[i]);
            exit(1);
        }
      else
        argv[j++] = argv[i];
    argc = j;

    if (argc != 2)
      { fprintf(stderr, "Usage: %s %s\n", Prog_Name, Usage);
        exit(1);
      }
  }

  gixpath = argv[1];
  dir = PathTo(gixpath);
  root = Root(gixpath, ".gix");

  // Read stub to get number of threads/parts
  { char *full = Malloc(strlen(dir) + strlen(root) + 20, "Path alloc");
    sprintf(full, "%s/%s.gix", dir, root);
    f = open(full, O_RDONLY);
    if (f < 0)
      { fprintf(stderr, "%s: Cannot open %s\n", Prog_Name, full);
        free(full);
        free(dir);
        free(root);
        exit(1);
      }
    int kmer;
    read(f, &kmer, sizeof(int));
    read(f, &nthreads, sizeof(int));
    close(f);
    free(full);
  }

  printf("%s %d ktab parts with zstd level %d, frame size %u\n",
         decompress ? "Decompressing" : "Compressing", nthreads, level, frameSize);

  // Process each ktab part
  for (p = 1; p <= nthreads; p++)
    { char inpath[1024], outpath[1024];

      if (decompress)
        { sprintf(inpath, "%s/.%s.ktab.%d.zst", dir, root, p);
          sprintf(outpath, "%s/.%s.ktab.%d", dir, root, p);
          printf("  %s -> %s\n", inpath, outpath);
          if (decompress_ktab(inpath, outpath) != 0)
            { fprintf(stderr, "Failed to decompress %s\n", inpath);
              exit(1);
            }
        }
      else
        { sprintf(inpath, "%s/.%s.ktab.%d", dir, root, p);
          sprintf(outpath, "%s/.%s.ktab.%d.zst", dir, root, p);
          printf("  %s -> %s\n", inpath, outpath);
          if (compress_ktab(inpath, outpath, level, frameSize) != 0)
            { fprintf(stderr, "Failed to compress %s\n", inpath);
              exit(1);
            }
        }
    }

  // Show compression ratio
  if (!decompress)
    { int64_t orig_total = 0, comp_total = 0;
      for (p = 1; p <= nthreads; p++)
        { char path[1024];
          struct stat st;
          sprintf(path, "%s/.%s.ktab.%d", dir, root, p);
          if (stat(path, &st) == 0)
            orig_total += st.st_size;
          sprintf(path, "%s/.%s.ktab.%d.zst", dir, root, p);
          if (stat(path, &st) == 0)
            comp_total += st.st_size;
        }
      printf("Original: %.2f MB, Compressed: %.2f MB (%.2fx)\n",
             orig_total / 1e6, comp_total / 1e6, (double)orig_total / comp_total);
    }

  free(dir);
  free(root);
  return 0;
}
