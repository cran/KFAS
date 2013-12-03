void F77_NAME(ldl)(double *a, int *n, double *tol, int *info);
void F77_NAME(ldlssm)(double *yt, int *ydimt, int *yobs, int *timevar, double *zt, int *p, int *m, int *n, double *ichols,int *nh,int *hchol,int *dim,int *info, int *hobs,double *tol);
void F77_NAME(signaltheta)(int *tvz, double *zt, double *ahat, double *vt, int *p, int *n, int *m, double *theta, double *thetavar,int *d,int *states,int *m2);
void F77_NAME(approx)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *ht, double *qt, double *a1, double *p1,double *p1inf, int *p,int *n,int *m,int *r, double *theta, double *u, double *ytilde, int *dist,int *maxiter,double *tol,int *rankp,double *convtol);
void F77_NAME(simgaussian)(int *ymiss,int *timevar, double *yt, double *zt, double *ht, double *tt, double *rtv, double *qt, double *a1, double *p1, double *p1inf, int *nnd, int *nsim, double *epsplus, double *etaplus, double *aplus1, int *p, int *n, int *m, int *r, int *info,int *rankp, double *tol, int *nd, int *ndl, double *sim, double *c, int *simwhat, int *simdim, int *antithetics);
void F77_NAME(gsmoothall)(int *ymiss, int *timevar, double *zt, double *ht,double *tt, double *rtv, double *qt, int *p, int *n, int *m, int *r, int *d,int *j, double *at, double *pt, double *vt, double *ft, double *kt, double *rt, double *rt0, double *rt1, double *nt, double *nt0, double *nt1, double *nt2, double *pinf, double *kinf,double *finf,  double *tol,double *ahat, double *vvt,double *epshat,double *epshatvar, double *etahat,double *etahatvar,double *thetahat,double *thetahatvar, int *ldlsignal,double *zorig, int *zorigtv, int *aug,int *state,int *dist,int *signal);
void F77_NAME(ngsmooth)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *qt, double *a1, double *p1,double *p1inf, double *u, double *theta,int *dist, int *p,int *n, int *m, int *r, int *rankp, int *nnd,int *nsim,double *epsplus,double *etaplus,double *aplus1,double *c,double *tol,int *info, int *maxiter,double *convtol,int *nd,int *ndl,double *alphahat,double *alphavar,double *thetahat,double *thetavar,double *yhat,double *yvar,int *smootha, int smooths, int smoothy);
void F77_NAME(kfilter)(double *yt, int *ymiss, int *timevar, double *zt, double *ht,double *tt, double *rt, double *qt, double *a1, double *p1, double *p1inf, int *p,int *n,int *m,int *r,int *d,int *j,double *at, double *pt, double *vt, double *ft,double *kt, double *pinf, double *finf, double *kinf, double *lik, double *tol,int *rankp, double *theta, double *thetavar, int *filtersignal);
void F77_NAME(glogliku)(double *yt, int *ymiss, int *timevar, double *zt, double *ht,double *tt, double *rt, double *qt, double *a1, double *p1, double *p1inf,int *m, int *r, int *n, double *lik, double *tol,int *rankp);
void F77_NAME(gloglik)(double *yt, int *ymiss, int *timevar, double *zt, double *ht, double *tt, double *rt, double *qt, double *a1, double *p1, double *p1inf,int *p, int *m, int *r, int *n, double *lik, double *tol,int *rankp);
void F77_NAME(ngloglik)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *qt, double *a1, double *p1,double *p1inf, int *p,int *m,int *r, int *n, double *lik, double *theta, double *u, int *dist,int *maxiter,int *rankp,double *convtol, int *nnd,int *nsim,double *epsplus,double *etaplus,double *aplus1,double *c,double *tol,int *info,int *antit,int *sim,int *nsim2,int *nd,int *ndl);
void F77_NAME(isample)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *qt, double *a1, double *p1,double *p1inf, double *u, int *dist, int *p, int *n, int *m, int *r, double *theta, int *maxiter,int *rankp,double *convtol, int *nnd,int *nsim, double *epsplus,double *etaplus,double *aplus1,double *c,double *tol,int *info,int *antithetics,double *w,double *sim,int *nd,int *ndl, int* simwhat, int* simdim);
void F77_NAME(zalpha)(int *timevar, double *zt, int *alpha, double *theta, int *p, int *m, int *n, int *nsim, int *m2, int *states);
void F77_NAME(varmeanw)(double *x,double *w,int *m,int *n,int *k,double *meanx,double *varx,int *var);
void F77_NAME(artransform)(double *u, double *phi, int *p);
void F77_NAME(simfilter)(int *ymiss,int *timevar, double *yt, double *zt, double *ht, double *tt, double *rtv, double *qt, double *a1, double *p1, double *p1inf, int *nnd, int *nsim, double *epsplus, double *etaplus, double *aplus1, int *p, int *n, int *m, int *r, int *info,int *rankp, double *tol, int *nd, int *ndl, double *sim, double *c, int *simwhat, int *simdim, int *antithetics);
void F77_NAME(ngfilter)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *qt, double *a1, double *p1,double *p1inf, double *u, double *theta,int *dist, int *p,int *n, int *m, int *r, int *rankp, int *nnd,int *nsim,double *epsplus,double *etaplus,double *aplus1,double *c,double *tol,int *info, int *maxiter,double *convtol,int *nd,int *ndl,double *alphahat,double *alphavar,double *thetahat,double *thetavar,double *yhat,double *yvar,int *smootha, int smooths, int smoothy);
void F77_NAME(isamplefilter)(double *yt, int *ymiss, int *timevar, double *zt, double *tt, double *rtv, double *qt, double *a1, double *p1,double *p1inf, double *u, int *dist, int *p, int *n, int *m, int *r, double *theta, int *maxiter,int *rankp,double *convtol, int *nnd,int *nsim, double *epsplus,double *etaplus,double *aplus1,double *c,double *tol,int *info,int *antithetics,double *w,double *sim,int *nd,int *ndl, int* simwhat, int* simdim);

