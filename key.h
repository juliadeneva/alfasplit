
#ifndef CHARSTAR
#include "y.tab.h" /* pull in the defs for the "type" enum */
#endif

struct HEADERKEY {
  void *next;
  char *name;
  int offset;
  int type;
  int len;  /* length of data element */
  int alen; /* array length */
};

struct HEADERVAL {
  void *value;
  struct HEADERKEY *key;
};

struct HEADERP {
  struct HEADERKEY *head;
  struct HEADERKEY *tail;
  char *buf;       /* ascii C header declaration */
  int offset;      /* len of buf ( offset to start of generic head in file */
  void *header;    /* pointer to instance of generic header */
  int headlen;     /* len of generic header */
  int fd;          /* file des */
  int yacc_offset; /* last returned by head_input */
};

struct HEADERP *head_parse(char *, int *);
extern struct HEADERKEY headerkey[];
