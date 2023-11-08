# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
@article{delaRosaTrevin2013,
title = "Xmipp 3.0: An improved software suite for image processing in electron microscopy ",
journal = "JSB",
volume = "184",
number = "2",
pages = "321 - 328",
year = "2013",
issn = "1047-8477",
doi = "http://dx.doi.org/10.1016/j.jsb.2013.09.015",
url = "http://www.sciencedirect.com/science/article/pii/S1047847713002566",
author = "de la Rosa-Trevín, J.M.  and Oton, J. and R. Marabini and A. Zaldívar and J. Vargas and J.M. Carazo and Sorzano, C.O.S.",
keywords = "Electron microscopy, Single particles analysis, Image processing, Software package "
}

@article{Jimenez2022,
title = {ScipionTomo: Towards cryo-electron tomography software integration, reproducibility, and validation},
journal = {Journal of Structural Biology},
volume = {214},
number = {3},
pages = {107872},
year = {2022},
issn = {1047-8477},
doi = {https://doi.org/10.1016/j.jsb.2022.107872},
url = {https://www.sciencedirect.com/science/article/pii/S1047847722000429},
author = {J. Jiménez de la Morena and P. Conesa and Y.C. Fonseca and F.P. {de Isidro-Gómez} and D. Herreros and E. Fernández-Giménez and D. Strelak and E. Moebel and T.O. Buchholz and F. Jug and A. Martinez-Sanchez and M. Harastani and S. Jonic and J.J. Conesa and A. Cuervo and P. Losana and I. Sánchez and M. Iceta and L. {del Cano} and M. Gragera and R. Melero and G. Sharov and D. Castaño-Díez and A. Koster and J.G. Piccirillo and J.L. Vilas and J. Otón and R. Marabini and C.O.S. Sorzano and J.M. Carazo},
abstract = {Image processing in cryogenic electron tomography (cryoET) is currently at a similar state as Single Particle Analysis (SPA) in cryogenic electron microscopy (cryoEM) was a few years ago. Its data processing workflows are far from being well defined and the user experience is still not smooth. Moreover, file formats of different software packages and their associated metadata are not standardized, mainly since different packages are developed by different groups, focusing on different steps of the data processing pipeline. The Scipion framework, originally developed for SPA (de la Rosa-Trevín et al., 2016), has a generic python workflow engine that gives it the versatility to be extended to other fields, as demonstrated for model building (Martínez et al., 2020). In this article, we provide an extension of Scipion based on a set of tomography plugins (referred to as ScipionTomo hereafter), with a similar purpose: to allow users to be focused on the data processing and analysis instead of having to deal with multiple software installation issues and the inconvenience of switching from one to another, converting metadata files, managing possible incompatibilities, scripting (writing a simple program in a language that the computer must convert to machine language each time the program is run), etcetera. Additionally, having all the software available in an integrated platform allows comparing the results of different algorithms trying to solve the same problem. In this way, the commonalities and differences between estimated parameters shed light on which results can be more trusted than others. ScipionTomo is developed by a collaborative multidisciplinary team composed of Scipion team engineers, structural biologists, and in some cases, the developers whose software packages have been integrated. It is open to anyone in the field willing to contribute to this project. The result is a framework extension that combines the acquired knowledge of Scipion developers in close collaboration with third-party developers, and the on-demand design of functionalities requested by beta testers applying this solution to actual biological problems.}
}

"""
