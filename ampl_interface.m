%AMPL_INTERFACE Interface to the AMPL solver library.
classdef ampl_interface < handle

   properties (SetAccess = private, Hidden = true)
      oH        % handle to the underlying C++ class instance
      interface % handle to dense or sparse interface
   end

   properties
      x0, bl, bu, v, cl, cu, nlc
      sigma     % scale of the Lagrangian
   end

   methods
      function self = ampl_interface(model, sparse)
         % Constructor   Create the ampl model instance.

         if nargin == 2 && sparse
            self.interface = @ampl.spam_interface_mex;
         else
            self.interface = @ampl.ampl_interface_mex;
         end

         [self.oH,...
            self.x0, self.bl, self.bu,...
            self.v, self.cl, self.cu, self.nlc] = ...
            self.interface('new', model);

         % Set the default sign of the Lagrangian.
         self.sigma = -1;
         self.lagscale(self.sigma);
      end

      function delete(self)
         % Destructor   Destroy the ampl model instance.
         if isempty(self.oH)
            % No need to delete if the constructor failed.
            return
         end
         self.interface('delete', self.oH);
      end

      function write_sol(self, msg, x, y)
         if issparse(x), x = full(x); end
         if issparse(y), y = full(y); end
         self.interface('write_sol', self.oH, msg, x, y);
      end

      function f = obj(self, x)
         if issparse(x), x = full(x); end
         f = self.interface('obj', self.oH, x);
      end

      function g = grad(self, x)
         if issparse(x), x = full(x); end
         g = self.interface('grad', self.oH, x);
      end

      function H = hessobj(self, x)  %#ok<INUSD>
         H = self.interface('hessobj', self.oH);
      end

      function c = con(self, x)
         if issparse(x), x = full(x); end
         c = self.interface('con', self.oH, x);
      end

      function J = jac(self, x)
         if issparse(x), x = full(x); end
         J = self.interface('jac', self.oH, x);
      end

      function HL = hesslag(self, y)
         if issparse(y), y = full(y); end
         HL = self.interface('hesslag', self.oH, y);
      end

      function Hv = hesslagprod(self, y, v)
         if issparse(y), y = full(y); end
         if issparse(v), v = full(v); end
         Hv = self.interface('hesslagprod', self.oH, y, v);
      end

      function HC = hesscon(self, y)
         if self.sigma < 0
            y = -y;
         end
         if issparse(y), y = full(y); end
         HC = self.interface('hesscon', self.oH, y);
      end

      function Hv = hessconprod(self, y, v)
         if self.sigma < 0
            y = -y;
         end
         if issparse(y), y = full(y); end
         Hv = self.interface('hessconprod', self.oH, y, v);
      end

      function gHiv = ghivprod(self, x, g, v)
         %GHIVPROD  Return m-vector of dot products <g,Hi*V>.
         if issparse(x), x = full(x); end
         if issparse(g), g = full(g); end
         if issparse(v), v = full(v); end
         gHiv = self.interface('ghivprod', self.oH, x, g, v);
      end

      function lagscale(self, scale)
         %LAGSCALE Set the scale of the Lagrangian.
         %
         % LAGSCALE(sigma) sets the scale of the Lagrangian as
         % L(x,y) = f(x) + sigma <c, y>.
         % By default, sigma = -1 when the object is created.
         self.interface('lagscale', self.oH, scale);
      end
   end
end
