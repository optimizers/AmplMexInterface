#include "class_handle.hpp"
#undef printf
#include "asl/aslinterface.h"
#include "mex.h"

// The class that we are interfacing to
class ASLMex {
public:
    ASL *asl;
    mwSize nvar, ncon, nnzj, nnzh;
private:
};


static double*
getDense(const mxArray *mp, const char *who, mwSize m)
{
    char msgbuf[256];
    mwSize m1, n1;

    if (mxIsSparse(mp)) {
        sprintf(msgbuf,"Expected %s to be a dense matrix",who);
        mexErrMsgTxt(msgbuf);
    }
    m1 = mxGetM(mp);
    n1 = mxGetN(mp);
    if (m1 != m || (n1 != 1 && m1)) {
        sprintf(msgbuf,
                "Expected %s to be %d x 1 rather than %d x %d\n",
                who, m, m1, n1);
        mexErrMsgTxt(msgbuf);
    }
    return mxGetPr(mp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    ASLMex* AC;
    ASL *asl;
    char *buf1, buf[512];

    fint err = 0;
    mwSize i, j, nvar, ncon, nlcon, nnzj, nnzh;

    // Get the command string.
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // -----------------------------------------------------------------
    // Command: New
    // -----------------------------------------------------------------
    if (!strcmp("new", cmd)) {

        if (mxGetString(prhs[1], buf1 = buf, sizeof(buf)))
            mexErrMsgTxt("Expected 'stub' as argument\n");

        // Return a handle to a new C++ instance
        AC = new ASLMex;

        asl = AC->asl = asl_init(buf);
        AC->nvar = nvar = asl_nvar(asl);
        AC->ncon = ncon = asl_ncon(asl);
        AC->nnzj = nnzj = asl_nnzj(asl);
        AC->nnzh = nnzh = asl_nnzh(asl);

        nlcon = asl_nlc(asl);
        plhs[7] = mxCreateDoubleScalar(nlcon);

        // Left-hand sides
        plhs[0] = convertPtr2Mat<ASLMex>(AC);
        double *x0p = mxGetPr(plhs[1] = mxCreateDoubleMatrix(nvar, 1, mxREAL));
        double *blp = mxGetPr(plhs[2] = mxCreateDoubleMatrix(nvar, 1, mxREAL));
        double *bup = mxGetPr(plhs[3] = mxCreateDoubleMatrix(nvar, 1, mxREAL));
        double *y0p = mxGetPr(plhs[4] = mxCreateDoubleMatrix(ncon, 1, mxREAL));
        double *clp = mxGetPr(plhs[5] = mxCreateDoubleMatrix(ncon, 1, mxREAL));
        double *cup = mxGetPr(plhs[6] = mxCreateDoubleMatrix(ncon, 1, mxREAL));

        // Obtain pointers to arrays inside ASL data structure
        double *x0 = asl_x0(asl);
        double *bl = asl_lvar(asl);
        double *bu = asl_uvar(asl);
        double *y0 = asl_y0(asl);
        double *cl = asl_lcon(asl);
        double *cu = asl_ucon(asl);

        // Copy arrays over to Matlab
        for (i = 0; i < AC->nvar; i++) {
          x0p[i] = x0[i];
          blp[i] = bl[i];
          bup[i] = bu[i];
        }
        for (i = 0; i < AC->ncon; i++) {
          y0p[i] = y0[i];
          clp[i] = cl[i];
          cup[i] = cu[i];
        }

        return;
    }

    // -----------------------------------------------------------------
    // Retrieve C++ object, and unpack it.
    // -----------------------------------------------------------------
    AC = convertMat2Ptr<ASLMex>(prhs[1]);
    asl = AC->asl;
    nvar = AC->nvar;
    ncon = AC->ncon;
    nnzj = AC->nnzj;
    nnzh = AC->nnzh;

    // -----------------------------------------------------------------
    // -----------------------------------------------------------------
    // Call the various class methods
    // -----------------------------------------------------------------
    // -----------------------------------------------------------------

    // -----------------------------------------------------------------
    // Command: Delete
    // -----------------------------------------------------------------
    if (!strcmp("delete", cmd)) {

        // Destroy the C++ object
        ASL_free(&(AC->asl));
        destroyObject<ASLMex>(prhs[1]);

        // Warn if other arguments were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // -----------------------------------------------------------------
    // Command: objective
    // -----------------------------------------------------------------
    if (!strcmp("obj", cmd)) {

        double *x = getDense(prhs[2], "x", nvar);
        double *f = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL));

        *f = asl_obj(asl, x, &err);
        if (err) mexErrMsgTxt("Trouble evaluating objective value\n");

        return;
    }

    // -----------------------------------------------------------------
    // Command: gradient
    // -----------------------------------------------------------------
    if (!strcmp("grad", cmd)) {

        double *x = getDense(prhs[2], "x", nvar);
        double *g = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nvar, 1, mxREAL));

        asl_grad(asl, x, g, &err);
        if (err) mexErrMsgTxt("Trouble evaluating objective gradient\n");

        return;
    }

    // -----------------------------------------------------------------
    // Command: constraint
    // -----------------------------------------------------------------
    if (!strcmp("con", cmd)) {

        double *x = getDense(prhs[2], "x", nvar);
        double *c = mxGetPr(plhs[0] = mxCreateDoubleMatrix(ncon, 1, mxREAL));

        asl_cons(asl, x, c, &err);
        if (err) mexErrMsgTxt("Trouble evaluating constraints\n");

        return;
    }

    // -----------------------------------------------------------------
    // Command: jacobian in coordinate format
    // -----------------------------------------------------------------
    if (!strcmp("jac_coord", cmd)) {

        double *x = getDense(prhs[2], "x", nvar);
        double *rows = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nnzj, 1, mxREAL));
        double *cols = mxGetPr(plhs[1] = mxCreateDoubleMatrix(nnzj, 1, mxREAL));
        double *vals = mxGetPr(plhs[2] = mxCreateDoubleMatrix(nnzj, 1, mxREAL));

        // is there a more elegant way?!
        int64_t *irows = (int64_t*)malloc(nnzj * sizeof(int64_t));
        int64_t *jcols = (int64_t*)malloc(nnzj * sizeof(int64_t));

        asl_jac(asl, x, irows, jcols, vals, &err);
        if (err) mexErrMsgTxt("Trouble evaluating constraints Jacobian\n");

        for (i = 0; i < nnzj; i++) {
          rows[i] = irows[i] + 1;
          cols[i] = jcols[i] + 1;  // 1-based indexing
        }
        free(irows);
        free(jcols);

        return;
    }

    // -----------------------------------------------------------------
    // Command: Hessian of Lagrangian
    // -----------------------------------------------------------------
    if (!strcmp("hesslag", cmd)) {

        double *y = getDense(prhs[2], "y", ncon);
        double ow = mxGetScalar(prhs[3]);

        double *rows = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nnzh, 1, mxREAL));
        double *cols = mxGetPr(plhs[1] = mxCreateDoubleMatrix(nnzh, 1, mxREAL));
        double *vals = mxGetPr(plhs[2] = mxCreateDoubleMatrix(nnzh, 1, mxREAL));

        int64_t *irows = (int64_t*)malloc(nnzh * sizeof(int64_t));
        int64_t *jcols = (int64_t*)malloc(nnzh * sizeof(int64_t));

        ow = asl->i.objtype_[0] ? -ow : ow;  // Objective weight.
        asl_hess(asl, y, ow, irows, jcols, vals);

        for (i = 0; i < nnzh; i++) {
          rows[i] = irows[i] + 1;
          cols[i] = jcols[i] + 1;
        }
        free(irows);
        free(jcols);

        return;
    }

    // -----------------------------------------------------------------
    // Command: lagscale
    // -----------------------------------------------------------------
    if (!strcmp("lagscale", cmd)) {

        // Specify that the sign of the Lagrangian is as follows:
        // L(x,y) = H(x) - sum_i y_i H_i(x).
        double sigma = mxGetScalar(prhs[2]);

        asl_lagscale(asl, sigma, &err);
        if (err) mexErrMsgTxt("Failed to set sign of Lagrangian.");

        return;
    }

    // -----------------------------------------------------------------
    // Command: ghivprod.  Vector of dot products <g, Hi*v>, where Hi
    // are Hessians of the constraints.
    // -----------------------------------------------------------------
    if (!strcmp("ghivprod", cmd)) {

        double *x = getDense(prhs[2], "x", nvar);
        double *g = getDense(prhs[3], "g", nvar);
        double *v = getDense(prhs[4], "v", nvar);
        // double *hv = (double*)mxMalloc(n*sizeof(double));
        double *gHiv = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nlc, 1, mxREAL));

        asl_ghjvprod(asl, g, v, gHiv);

        return;
    }

    // -----------------------------------------------------------------
    // Command: hesslagprod.  (H - sum_i H_i y_i) v
    // -----------------------------------------------------------------
    if (!strcmp("hesslagprod", cmd)) {

        double *y = getDense(prhs[2], "y", ncon);
        double *v = getDense(prhs[3], "v", nvar);
        double ow = mxGetScalar(prhs[4]);
        ow = asl->i.objtype_[0] ? -ow : ow;  // Objective weight.
        double *hv = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nvar, 1, mxREAL));

        asl_hprod(asl, y, v, hv, ow);

        return;
    }

    // -----------------------------------------------------------------
    // Command: write_sol.  Write solution to file.
    // -----------------------------------------------------------------
    if (!strcmp("write_sol", cmd)) {

        if (mxGetString(prhs[2], buf, sizeof(buf)))
            mexErrMsgTxt("Error while retrieving message.");
        double *x = getDense(prhs[3], "x", nvar);
        double *y = getDense(prhs[4], "y", ncon);

        asl_write_sol(asl, buf, x, y);

        return;
    }

    // -----------------------------------------------------------------
    // Got here, so command not recognized
    // -----------------------------------------------------------------
    mexErrMsgTxt("Command not recognized.");
}
