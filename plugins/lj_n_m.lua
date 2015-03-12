-- Copyright (c) 2011-2015, Florent Hedin, Markus Meuwly, and the University of Basel
-- All rights reserved.
-- The 3-clause BSD license is applied to this software.
-- see LICENSE.txt
-- 
-- lj_n_m.lua
--
-- Provides energy and gradient functions
-- for a modified n-m Lennard Jones potential
--
-- Formula is :
-- V = C * epsilon * [ (sigma/r)^n - (sigma/r)^m ]
-- Where n and m are integers and n>m ; C is a constant depending on both n and m
--
-- See for example : http://www.sklogwiki.org/SklogWiki/index.php/Lennard-Jones_model#n-m_Lennard-Jones_potential
--
-- at least 4 times faster when using luajit (see http://luajit.org/ )
--

-- parameters to tune for adapting the potential
-- this set is the same than the standard 12,6 potential
local c = 4.0
local n = 12.0
local m = 6.0

-- an example of parameters for the 9,6 parameters
-- see http://www.sklogwiki.org/SklogWiki/index.php/9-6_Lennard-Jones_potential
-- local c = 6.75
-- local n = 9.0
-- local m = 6.0

-- an example of parameters for the 9,3 parameters
-- see http://www.sklogwiki.org/SklogWiki/index.php/9-3_Lennard-Jones_potential
-- local c = 3.0*math.sqrt(3.0)/2.0
-- local n = 9.0
-- local m = 3.0

-- This function estimates the LJ n-m potential between a pair of atoms a and b
function lj_v_n_m_pair(xa, ya, za, xb, yb, zb, eps_a, eps_b, sig_a, sig_b)
    
    -- Lorentz-Berthelot rules
    local eps = math.sqrt(eps_a*eps_b)
    local sig = 0.5*(sig_a + sig_b)
    
    -- distance calculation
    local dx = xb - xa
    local dy = yb - ya
    local dz = zb - za
    local r = math.sqrt(dx*dx + dy*dy + dz*dz)
    
    -- evaluate pair potential
    local pot = c*eps*( math.pow(sig/r,n) - math.pow(sig/r,m) )
    
    return pot
    
end

---------------------------------------------------------------------------------

-- This function estimates the LJ n-m gradient contribution of one atom on another one
function lj_dv_n_m_pair(xa, ya, za, xb, yb, zb, eps_a, eps_b, sig_a, sig_b)
    
    -- Lorentz-Berthelot rules
    local eps = math.sqrt(eps_a*eps_b)
    local sig = 0.5*(sig_a + sig_b)
    
    -- distance calculation
    local dx = xa - xb
    local dy = ya - yb
    local dz = za - zb
    local r2 = dx*dx + dy*dy + dz*dz
    local r = math.sqrt(r2)
    
    -- evaluate gradient
    local dv = -c*eps*( n*math.pow(sig/r,n) - m*math.pow(sig/r,m) )/r2
    local fx = dv*dx
    local fy = dv*dy
    local fz = dv*dz
    
    -- note that we return fz first and fx last, as lua is stack based.
    -- By doing that the c code will take first fx as it is on top of the stack
    return fz,fy,fx
    
end



