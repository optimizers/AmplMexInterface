classdef tests < matlab.unittest.TestCase

   properties (TestParameter)
      hsproblem = {...
         'hs006', 'hs007', 'hs008', 'hs010', 'hs011', 'hs012', 'hs013',...
         'hs014', 'hs015', 'hs016', 'hs017', 'hs018', 'hs019', 'hs020',...
         'hs022', 'hs023', 'hs026', 'hs027', 'hs029', 'hs030', 'hs031',...
         'hs032', 'hs033', 'hs034', 'hs039', 'hs040', 'hs042', 'hs043',...
         'hs046', 'hs047', 'hs056', 'hs057', 'hs059', 'hs060', 'hs061',...
         'hs063', 'hs064', 'hs065', 'hs066', 'hs070', 'hs071', 'hs072',...
         'hs073', 'hs074', 'hs075', 'hs077', 'hs078', 'hs079', 'hs080',...
         'hs081', 'hs083', 'hs084', 'hs088', 'hs089', 'hs090', 'hs091',...
         'hs092', 'hs093', 'hs095', 'hs096', 'hs097', 'hs098', 'hs099',...
         'hs100', 'hs100lnp', 'hs100mod', 'hs102', 'hs103', 'hs104',...
         'hs106', 'hs107', 'hs108', 'hs109', 'hs111', 'hs111lnp', 'hs113',...
         'hs114', 'hs116', 'hs117', 'hs99exp'} ;
   end
   
   methods (Test)
      
      function test_ghiv_prod(testCase, hsproblem) 
         import matlab.unittest.constraints.RelativeTolerance;
         import matlab.unittest.constraints.IsEqualTo;

         rng('default')
         fname = ['/Users/mpf/local/TestProblems/compiled/' hsproblem];
         p = amplmodel(fname);
         
         % Manually
         nlc = sum(~p.linear);
         n = p.n; m = p.m;
         gHiv = zeros(nlc,1);
         g = randn(n,1);
         v = randn(n,1);
         y = zeros(m,1);
         for i=1:nlc
             y(i) = 1;
             Hi = p.hesscon(p.x0, y);
             y(i) = 0;
             gHiv(i) = g'*Hi*v;
         end
         
         % Using ghiv_prod
         gHiv2 = p.ghivprod(p.x0, g, v);
         
         testCase.verifyThat(gHiv, IsEqualTo(gHiv2, ...
            'Within', RelativeTolerance(sqrt(eps))));
         
      end
      
   end
   
end % classdef

