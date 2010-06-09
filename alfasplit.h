/* alfasplit.c */
int main(int argc, char *argv[]);
int new_filename(char *old, struct HEADERP *h, char *name[]);
int find_index(char *s);
int find_wappno(char *s);
int obs2mjd(char *dt);
void fix_positions(double ra, double dec, struct HEADERP *h, int beam);
double deg2dms(double t);
double dms2deg(double in);
int ast_seconds(char *s);
double rotx(double offx, double offy, double ang);
double roty(double offx, double offy, double ang);
