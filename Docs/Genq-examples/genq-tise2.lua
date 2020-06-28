-- Generalized coordinates for TiSe2


genq = { }

genq.nqs = 2
genq.natoms = 12
genq.q = { }
genq.xa = { }
genq.species = { }

function coords_from_qs ( q )
       -- maybe 'assert' that size(q)= genq.nqs
       
	 local delta1 = q[1]
	 local delta2 = q[2]

	 local xTi = 1/2
	 local xSe = 1/6
	 local ySe = 5/6
	 local zSe = 0.030402381   -- hard-wired value
	 
       local Ti_3e_1 = { xTi + delta1 , 0.0, 0.0 }
       local Ti_3e_2 = { 0 ,  xTi + delta1 , 0.0 }
       local Ti_3e_3 = { -xTi - delta1 , -xTi - delta1, 0.0 }

       local Se_2d_1 = { 1/3, 2/3, zSe }
       local Se_2d_2 = { 2/3, 1/3, -zSe }

       local Se_6g_1 = {  xSe - delta2 ,  ySe - delta2,  zSe }
       local Se_6g_2 = { -ySe + delta2 ,  xSe - ySe   ,  zSe }
       local Se_6g_3 = { -xSe + ySe    , -xSe + delta2,  zSe }
       local Se_6g_4 = {  ySe - delta2 ,  xSe - delta2, -zSe }
       local Se_6g_5 = {  xSe - ySe    , -ySe + delta2, -zSe }
       local Se_6g_6 = { -xSe + delta2 ,  xSe + ySe   , -zSe }

       local xa  = { Ti_3e_1, Ti_3e_2,
                     Se_2d_1, Se_2d_2,
		     Se_6g_1, Se_6g_2, Se_6g_3,
		     Se_6g_4, Se_6g_5, Se_6g_6 }
		     
       local isa = { 1,  1,
                     2,  2,
		     2,  2, 2, 2, 2, 2}

       return xa, isa
       -- maybe assert that size of tables  == genq.natoms
end

