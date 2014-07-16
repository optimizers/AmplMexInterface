#include "class_handle.hpp"
#undef printf
#include "asl/asl_pfgh.h"
#include "mex.h"

// The class that we are interfacing to
class dummy {
public:
    ASL *asl;
    real *J;
    fint n, nc, nz;
    size_t Jsize;
private:
};

extern "C" {
    double ddot_(
    size_t *n,
    double *dx,
    size_t *incx,
    double *dy,
    size_t *incy
);
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
    dummy* AC;
    ASL *asl;
    char *buf1, buf[512], msgbuf[256];
    int nerror = 0;
    mwSize i, n, nc, nz;

    // Get the command string.
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // -----------------------------------------------------------------
    // Command: New
    // -----------------------------------------------------------------
    if (!strcmp("new", cmd)) {

        FILE *nl;

        // Return a handle to a new C++ instance
        AC = new dummy;
        plhs[0] = convertPtr2Mat<dummy>(AC);

        // Allocate the ASL object.
        if (mxGetString(prhs[1], buf1 = buf, sizeof(buf)))
            mexErrMsgTxt("Expected 'stub' as argument\n");
        asl = AC->asl = ASL_alloc(ASL_read_pfgh);
        return_nofile = 1;
        if (!(nl = jac0dim(buf1,strlen(buf)))) {
            sprintf(msgbuf, "Can't open %.*s\n", sizeof(msgbuf)-20, buf);
            mexErrMsgTxt(msgbuf);
        }
        if (n_obj <= 0)
            printf("Warning: objective == 0\n");

        // Save the sizes into the ampl context. The sizes n, nc, and
        // nz are useful, so make local copies.
        n = AC->n = n_var;
        nc = AC->nc = n_con;
        nz = AC->nz = nzc;

        // Allocate space for the Jacobian. M1alloc is an ampl macro
        // that actually allocates space within the ASL object; no
        // need to deallocate this, since the ASL_free function
        // (later) will take care to delete this.
        AC->J = (real *)M1alloc(nz*sizeof(real));

        // Create the RHSs.
        X0    = mxGetPr(plhs[1] = mxCreateDoubleMatrix(n , 1, mxREAL));
        LUv   = mxGetPr(plhs[2] = mxCreateDoubleMatrix(n , 1, mxREAL));
        Uvx   = mxGetPr(plhs[3] = mxCreateDoubleMatrix(n , 1, mxREAL));
        pi0   = mxGetPr(plhs[4] = mxCreateDoubleMatrix(nc, 1, mxREAL));
        LUrhs = mxGetPr(plhs[5] = mxCreateDoubleMatrix(nc, 1, mxREAL));
        Urhsx = mxGetPr(plhs[6] = mxCreateDoubleMatrix(nc, 1, mxREAL));
        plhs[7] = mxCreateDoubleScalar(nlc);

        // pfgh_read reads and then closes the open file nl.
        pfgh_read(nl, ASL_findgroups);
        AC->Jsize = nc*n*sizeof(real);

        return;
    }

    // -----------------------------------------------------------------
    // Retrieve C++ object, and unpack it.
    // -----------------------------------------------------------------
    AC = convertMat2Ptr<dummy>(prhs[1]);
    asl = AC->asl;
    n = AC->n;
    nc = AC->nc;
    nz = AC->nz;

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
        destroyObject<dummy>(prhs[1]);

        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // -----------------------------------------------------------------
    // Command: objective
    // -----------------------------------------------------------------
    if (!strcmp("obj", cmd)) {
        double *x = getDense(prhs[2], "x", n);
        double *f = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL));
        *f = objval(0, x, &nerror);
        if (nerror)
            mexErrMsgTxt("Trouble evaluating f\n");
        return;
    }

    // -----------------------------------------------------------------
    // Command: gradient
    // -----------------------------------------------------------------
    if (!strcmp("grad", cmd)) {
        double *x = getDense(prhs[2], "x", n);
        double *g = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL));
        objgrd(0, x, g, &nerror);
        if (nerror)
            mexErrMsgTxt("Trouble evaluating g\n");
        return;
    }

    // -----------------------------------------------------------------
    // Command: constraint
    // -----------------------------------------------------------------
    if (!strcmp("con", cmd)) {

        double *x = getDense(prhs[2], "x", n);
        double *c = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nc, 1, mxREAL));
        conval(x, c, &nerror);
        if (nerror)
            mexErrMsgTxt("Trouble evaluating c\n");
        return;
    }

    // -----------------------------------------------------------------
    // Command: jacobian
    // -----------------------------------------------------------------
    if (!strcmp("jac", cmd)) {

        cgrad *cg, **cgp, **cgpe;
        double *x = getDense(prhs[2], "x", n);
        double *J1 = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nc, n, mxREAL));
        if (nc) {
            memset(J1, 0, AC->Jsize);
            jacval(x, AC->J, &nerror);
            if (nerror)
                mexErrMsgTxt("Trouble evaluating J\n");
            cgp = Cgrad;
            for(cgpe = cgp + nc; cgp < cgpe; J1++)
                for(cg = *cgp++; cg; cg = cg->next)
                    J1[nc*cg->varno] = AC->J[cg->goff];
        }
        return;
    }

    // -----------------------------------------------------------------
    // Command: Hessian of Lagrangian
    // -----------------------------------------------------------------
    if (!strcmp("hessobj", cmd)) {

        double *H = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL));
	fullhes(H, n, 0, 0, NULL);
        return;
    }

    // -----------------------------------------------------------------
    // Command: Hessian of Lagrangian
    // -----------------------------------------------------------------
    if (!strcmp("hesslag", cmd)) {

        double *y = getDense(prhs[2], "y", nc);
        double *H = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL));
	fullhes(H, n, 0, 0, y);
        return;
    }

    // -----------------------------------------------------------------
    // Command: Hessian of constraint, sum_i H_i y_i.
    // -----------------------------------------------------------------
    if (!strcmp("hesscon", cmd)) {

        double *y = getDense(prhs[2], "y", nc);
        double *H = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL));
	fullhes(H, n, -1, NULL, y);
        return;
    }

    // -----------------------------------------------------------------
    // Command: lagscale
    // -----------------------------------------------------------------
    if (!strcmp("lagscale", cmd)) {

        // Specify that the sign of the Lagrangian is as follows:
        // L(x,y) = H(x) - sum_i y_i H_i(x).
        double sigma = mxGetScalar(prhs[2]);
        lagscale(sigma, &nerror);
        if (nerror)
            mexErrMsgTxt("Failed to set sign of Lagrangian Hessian.");
        return;
    }

    // -----------------------------------------------------------------
    // Command: ghivprod.  Vector of dot products <g, Hi*v>, where Hi
    // are Hessians of the constraints.
    // -----------------------------------------------------------------
    if (!strcmp("ghivprod", cmd)) {

        double *x = getDense(prhs[2], "x", n);
        double *g = getDense(prhs[3], "g", n);
        double *v = getDense(prhs[4], "v", n);
        double *hv = (double*)mxMalloc(n*sizeof(double));
        double *gHiv = mxGetPr(plhs[0] = mxCreateDoubleMatrix(nlc, 1, mxREAL));
        size_t one = 1;
        size_t nn = n;
        xknown(x);
        for (i = 0; i < nlc; i++) {
            hvcompd(hv, v, i);
            gHiv[i] = ddot_(&nn, hv, &one, g, &one);
        }
        xunknown();
        mxFree(hv);
        return;
    }

    // -----------------------------------------------------------------
    // Got here, so command not recognized
    // -----------------------------------------------------------------
    mexErrMsgTxt("Command not recognized.");
}
