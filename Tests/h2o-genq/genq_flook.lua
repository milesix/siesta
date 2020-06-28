-- Generalized coordinates for water molecule


genq = { }

genq.nqs = 2
genq.natoms = 3
genq.q = { }
genq.xa = { }
genq.species = { }

function coords_from_qs ( q )
       -- maybe 'assert' that size(q)= genq.nqs
       
	 local doh = q[1]
	 local alpha = q[2]

	 local sin = math.sin
	 local cos = math.cos
	 
       local H1 = { doh , 0.0, 0.0 }
       local H2 = { doh*cos(alpha) , doh*sin(alpha), 0.0 }
       local O  = { 0.0 , 0.0, 0.0 }

       local xa  = { O, H1, H2 }
       local isa = { 1,  2,  2 }

       return xa, isa
       -- maybe assert that size of tables  == genq.natoms
end

