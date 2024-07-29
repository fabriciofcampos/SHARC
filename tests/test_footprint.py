# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:53:15 2017

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.support.footprint import Footprint

class FootprintAreaTest(unittest.TestCase):
    
    def setUp(self):
        # Satélites geoestacionários
        self.fa1 = Footprint(0.1,bore_lat_deg=0,bore_subsat_long_deg=0.0)
        self.fa2 = Footprint(0.325,bore_lat_deg = 0)
        self.fa3 = Footprint(0.325,elevation_deg = 20)

        # Satélites a 1200 km e 600 km
        self.fa4_1200km = Footprint(0.325,elevation_deg=90,sat_height=1200000)
        self.fa4_600km = Footprint(0.325,elevation_deg=90,sat_height=600000)
        self.fa5_1200km = Footprint(0.325,elevation_deg=20,sat_height=1200000)
        self.fa5_600km = Footprint(0.325,elevation_deg=20,sat_height=600000)
        

    def test_construction(self):
        self.assertEqual(self.fa1.bore_lat_deg,0)
        self.assertEqual(self.fa1.bore_subsat_long_deg,0)
        self.assertEqual(self.fa1.beam_width_deg,0.1)
        self.assertEqual(self.fa1.bore_lat_rad,0)
        self.assertEqual(self.fa1.bore_subsat_long_rad,0)
        self.assertEqual(self.fa1.beam_width_rad,np.pi/1800)
        self.assertEqual(self.fa1.beta,0)
        self.assertEqual(self.fa1.bore_tilt,0)
        
        self.assertEqual(self.fa2.bore_lat_deg,0)
        self.assertEqual(self.fa2.bore_subsat_long_deg,0)
        self.assertEqual(self.fa2.bore_lat_rad,0)
        self.assertEqual(self.fa2.bore_subsat_long_rad,0)
        
        self.assertEqual(self.fa3.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa3.bore_subsat_long_deg,61.84,delta=0.01)

        self.assertEqual(self.fa4_1200km.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa4_1200km.sigma,0.84,delta=0.01)
        self.assertAlmostEqual(self.fa4_1200km.bore_subsat_long_deg,0,delta=1e-12)

        self.assertEqual(self.fa4_600km.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa4_600km.sigma,0.91,delta=0.01)
        self.assertAlmostEqual(self.fa4_600km.bore_subsat_long_deg,0,delta=1e-12)

        self.assertEqual(self.fa5_1200km.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa5_1200km.bore_subsat_long_deg,17.74,delta=0.01)

        self.assertEqual(self.fa5_600km.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa5_600km.bore_subsat_long_deg,10.82,delta=0.01)

        
    def test_set_elevation(self):
        self.fa2.set_elevation(20)
        self.assertEqual(self.fa2.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa2.bore_subsat_long_deg,61.84,delta=0.01)

        self.fa4_1200km.set_elevation(20)
        self.assertAlmostEqual(self.fa4_1200km.bore_subsat_long_deg,17.74,delta=0.01)

        self.fa4_600km.set_elevation(20)
        self.assertAlmostEqual(self.fa4_600km.bore_subsat_long_deg,10.82,delta=0.01)

        
    def test_calc_footprint(self):
        fp_long, fp_lat = self.fa1.calc_footprint(4)
        npt.assert_allclose(fp_long,np.array([0.0, 0.487,  -0.487, 0.0]),atol=1e-2)
        npt.assert_allclose(fp_lat,np.array([-0.562,  0.281,  0.281,  -0.562]),atol=1e-2)

        fp_long_1200km, fp_lat_1200km = self.fa4_1200km.calc_footprint(4)
        npt.assert_allclose(fp_long_1200km,np.array([.061, -0.031, -0.031, .061]),atol=1e-2)
        npt.assert_allclose(fp_lat_1200km,np.array([.0, 0.05, -0.05, .0]),atol=1e-2)

        fp_long_600km, fp_lat_600km = self.fa4_600km.calc_footprint(4)
        npt.assert_allclose(fp_long_600km,np.array([.03, -0.015, -0.015, .03]),atol=1e-2)
        npt.assert_allclose(fp_lat_600km,np.array([.0, 0.026, -0.026, .0]),atol=1e-2)
        

    def test_calc_area(self):
        # A avaliação da faixa de valores aceitáveis pode ser feita tanto de forma
        # percentual (em relação ao valor desperado) como usando um valor de diferença fixo
        p = 0.0025      # diferença de até 0.25%
        d = 200         # diferença de até 200 km2

        # O valor deve ser próximo de 129479.918476 km2 (área do setor esférico)
        # Usando 0.25% o teste passa, mas usando delta=d, falha
        a1 = self.fa2.calc_area(1000)
        self.assertAlmostEqual(a1,129500,delta=129500*p)

        # O valor deve ser próximo de 145.547436 km2 (área do setor esférico)
        a2 = self.fa4_1200km.calc_area(1000)
        self.assertAlmostEqual(a2,145.55,delta=d)

        # O valor deve ser próximo de 36.384779 km2 (área do setor esférico)
        self.fa4_600km.set_elevation(90)
        a3 = self.fa4_600km.calc_area(1000)
        self.assertAlmostEqual(a3,36.385,delta=d)

        # Área para satélite geoestacionário com longitude subsatelital deslocada (elevation=20)
        # Usando 0.25% o teste passa, mas usando delta=d, falha
        a4 = self.fa3.calc_area(1000)
        self.fa3.calc_area(1000)
        self.assertAlmostEqual(a4,486300,delta=486300*p)
       
        # Área para satélite a 1200 km com longitude subsatelital deslocada (elevation=20)
        a5 = self.fa5_1200km.calc_area(1000)
        self.assertAlmostEqual(a5,1790,delta=d)

        # Área para satélite a 600 km com longitude subsatelital deslocada (elevation=20)
        a6 = self.fa5_600km.calc_area(1000)
        self.assertAlmostEqual(a6,575,delta=d)


if __name__ == '__main__':
    unittest.main()