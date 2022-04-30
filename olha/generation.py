#!/usr/bin/env python3

import olga
import sequence_generation
import tempfile
from olga.utils import calc_steady_state_dist
from olha.utils import gene_map, igor_genes_to_idx
import numpy as np
import random


class SequenceGeneration():
    """ Class that defines the null distribution of the sequences
    @Members:
    * genomic_data
    * generative_model
    * error_rate
    * distance
    * Vs/Js
    """

    def __init__(self, genomic_data, generative_model,
                 Vs=None, Js=None, error_rate=0., seed=None):
        """ Load the model, and generate sequence
        @ Arguments:
        * genomic_data, generative_model: olga model
        * Vs/Js: list considered V genes / J genes (default None:
        all V/J genes), there some layway on the gene name, for ex
        TRBV1 will use all TRBV1-* genes.
        * error_rate: per base pair error rate for the sequences
        """
        if seed is not None:
            random.seed(seed)

        self.error_rate = error_rate
        self.genomic_data = genomic_data
        self.generative_model = generative_model
        recombination_type = "VDJ" if hasattr(genomic_data, "genD") else "VJ"

        dctVs = igor_genes_to_idx(self.genomic_data, "V")
        dctJs = igor_genes_to_idx(self.genomic_data, "J")

        self.gen = None
        with tempfile.NamedTemporaryFile() as tmp:
            # write a temporary file in Igor format
            if recombination_type == "VDJ":
                if Vs is None:
                    self.iVs = range(self.generative_model.PV.shape[0])
                else:
                    self.iVs = [dctVs[vv]
                                for v in Vs for vv in gene_map(v, self.genomic_data)]
                if Js is None:
                    self.iJs = range(self.generative_model.PDJ.shape[1])
                else:
                    self.iJs = [dctJs[jj]
                                for j in Js for jj in gene_map(j, self.genomic_data)]
                if not self.write_restricted_recombination_model_VDJ(tmp.name):
                    Exception("Error during model creation")
                self.gen = sequence_generation.SequenceGenerationVDJ(tmp.name,
                                                                     random.randrange(
                                                                         0, 1000000000),
                                                                     False)
            elif recombination_type == "VJ":
                if Vs is None:
                    self.iVs = range(self.generative_model.PVJ.shape[0])
                else:
                    self.iVs = [dctVs[vv] for v in Vs
                                for vv in gene_map(v, self.genomic_data)]
                if Js is None:
                    self.iJs = range(self.generative_model.PVJ.shape[1])
                else:
                    self.iJs = [dctJs[jj] for j in Js
                                for jj in gene_map(j, self.genomic_data)]
                if not self.write_restricted_recombination_model_VJ(tmp.name):
                    Exception("Error during model creation")
                self.gen = sequence_generation.SequenceGenerationVJ(tmp.name,
                                                                    random.randrange(
                                                                        0, 1000000000),
                                                                    False)
            else:
                print('Unrecognized recombination type, should be "VDJ" or "VJ"')

        # make the inverse mapping for V/J
        if Vs is None:
            self.invV = {ii: v[0]
                         for ii, v in enumerate(self.genomic_data.genV)}
        else:
            self.invV = {ii: vname for ii, vname
                         in enumerate(vv for v in Vs
                                      for vv in gene_map(v, self.genomic_data)
                                      )}
        if Js is None:
            self.invJ = {ii: j[0]
                         for ii, j in enumerate(self.genomic_data.genJ)}
        else:
            self.invJ = {ii: jname for ii, jname
                         in enumerate(jj for j in Js
                                      for jj in gene_map(j, self.genomic_data)
                                      )}

    def gen_rnd_prod_CDR3(self):
        """ Use Olga to generate a valid sequence with the
        generation probabilities of the model inferred by
        IGoR. Olga does not add errors, so we add them by
        hand.
        @ Return:
        (nt, aa, V, J)
        """
        nt, aa, V, J = self.gen.generate(True)
        return nt, aa, self.invV[V], self.invJ[J]

    def gen_rnd_CDR3(self):
        """ Use Olga to generate a valid sequence with the
        generation probabilities of the model inferred by
        IGoR. Olga does not add errors, so we add them by
        hand.
        @ Return:
        (nt, aa, V, J)
        """
        nt, aa, V, J = self.gen.generate(False)
        return nt, aa, self.invV[V], self.invJ[J]

    def write_restricted_recombination_model_VDJ(self,
                                                 filename):
        """ Write a file containing all the data related with the
        recombination model.
        Return true if everything went well.
        """
        PV = self.generative_model.PV
        PDJ = self.generative_model.PDJ
        PinsVD = self.generative_model.PinsVD
        PinsDJ = self.generative_model.PinsDJ
        PdelV_given_V = self.generative_model.PdelV_given_V
        PdelJ_given_J = self.generative_model.PdelJ_given_J
        PdelDldelDr_given_D = self.generative_model.PdelDldelDr_given_D

        PV = PV[self.iVs]

        if(np.sum(PV) == 0.):
            return False

        PV /= np.sum(PV)
        PDJ = PDJ[:, self.iJs]

        if(np.sum(PDJ) == 0.):
            return False

        PDJ /= np.sum(PDJ)

        PdelV_given_V = self.generative_model.PdelV_given_V[:,
                                                            self.iVs]
        PdelJ_given_J = self.generative_model.PdelJ_given_J[:,
                                                            self.iJs]

        if self.generative_model.first_nt_bias_insVD is None:
            first_nt_bias_insVD = calc_steady_state_dist(
                self.generative_model.Rvd)
        else:
            first_nt_bias_insVD = self.generative_model.first_nt_bias_insVD

        if self.generative_model.first_nt_bias_insDJ is None:
            first_nt_bias_insDJ = calc_steady_state_dist(
                self.generative_model.Rdj)
        else:
            first_nt_bias_insDJ = self.generative_model.first_nt_bias_insDJ

        error_rate = self.error_rate

        # Write the file
        with open(filename, "w") as fw:
            # Start by writing out the V, D and J genes
            print("# V genes", file=fw)
            print(len(self.iVs), file=fw)
            for iV in self.iVs:
                print(self.genomic_data.cutV_genomic_CDR3_segs[iV], file=fw)

            print("# D genes", file=fw)
            print(len(self.genomic_data.cutD_genomic_CDR3_segs), file=fw)
            for D in self.genomic_data.cutD_genomic_CDR3_segs:
                print(D, file=fw)

            print("# J genes", file=fw)
            print(len(self.iJs), file=fw)
            for iJ in self.iJs:
                print(self.genomic_data.cutJ_genomic_CDR3_segs[iJ], file=fw)

            # Now write all the probabilities
            print("# Marginals", file=fw)
            print("# P(V)", file=fw)
            print(len(PV), file=fw)
            print(*PV, file=fw)

            print("# P(D, J)", file=fw)
            print(*PDJ.shape, file=fw)
            for x in PDJ:
                print(*x, file=fw, end=" ")
            print(file=fw)

            print("# P(delV | V)", file=fw)
            print(*PdelV_given_V.transpose().shape, file=fw)
            for x in PdelV_given_V.transpose():
                print(*x, file=fw)

            print("# P(delJ | J)", file=fw)
            print(*PdelJ_given_J.transpose().shape, file=fw)
            for x in PdelJ_given_J.transpose():
                print(*x, file=fw)

            print("# P(delDl, delDr | D) [P(delD5, delD3 | D)]", file=fw)
            print(*PdelDldelDr_given_D.shape, file=fw)
            for y in PdelDldelDr_given_D.transpose():
                x = y.transpose()
                for a in x:
                    print(*a, file=fw, end=" ")
                print(file=fw)

            print("# P(insVD)", file=fw)
            print(len(PinsVD), file=fw)
            print(*PinsVD, file=fw)

            print("# Markov VD", file=fw)
            for i in range(4):
                for j in range(4):
                    print(self.generative_model.Rvd[j, i], file=fw, end=" ")
            print(file=fw)

            print("# First nucleotide bias VD", file=fw)
            print(*first_nt_bias_insVD, file=fw)

            print("# P(insDJ)", file=fw)
            print(len(PinsDJ), file=fw)
            print(*PinsDJ, file=fw)

            print("# Markov DJ", file=fw)
            for i in range(4):
                for j in range(4):
                    print(self.generative_model.Rdj[j, i], file=fw, end=" ")
            print(file=fw)

            print("# First nucleotide bias DJ", file=fw)
            print(*first_nt_bias_insDJ, file=fw)

            print("# Error rate", file=fw)
            print(error_rate, file=fw)

            print("# Thymus selection parameter", file=fw)
            print(1., file=fw)

            print("# Conserved J residues", file=fw)
            print("F V W", file=fw)

        return True

    def write_restricted_recombination_model_VJ(self,
                                                filename):
        """ Write a file containing all the data related with the
        recombination model (VJ version).
        Return true if everything went well.
        """
        PVJ = self.generative_model.PVJ
        PinsVJ = self.generative_model.PinsVJ
        PdelV_given_V = self.generative_model.PdelV_given_V
        PdelJ_given_J = self.generative_model.PdelJ_given_J

        PVJ = PVJ[np.ix_(self.iVs, self.iJs)]

        if(np.sum(PVJ) == 0.):
            return False

        PVJ /= np.sum(PVJ)

        PdelV_given_V = self.generative_model.PdelV_given_V[:,
                                                            self.iVs]
        PdelJ_given_J = self.generative_model.PdelJ_given_J[:,
                                                            self.iJs]

        if self.generative_model.first_nt_bias_insVJ is None:
            first_nt_bias_insVJ = calc_steady_state_dist(
                self.generative_model.Rvj)
        else:
            first_nt_bias_insVJ = self.generative_model.first_nt_bias_insVJ

        error_rate = self.error_rate

        # Write the file
        with open(filename, "w") as fw:
            # Start by writing out the V, D and J genes
            print("# V genes", file=fw)
            print(len(self.iVs), file=fw)
            for iV in self.iVs:
                print(self.genomic_data.cutV_genomic_CDR3_segs[iV], file=fw)

            print("# J genes", file=fw)
            print(len(self.iJs), file=fw)
            for iJ in self.iJs:
                print(self.genomic_data.cutJ_genomic_CDR3_segs[iJ], file=fw)

            # Now write all the probabilities
            print("# P(V, J)", file=fw)
            print(*PVJ.shape, file=fw)
            for x in PVJ:
                print(*x, file=fw, end=" ")
            print(file=fw)

            print("# P(delV | V)", file=fw)
            print(*PdelV_given_V.transpose().shape, file=fw)
            for x in PdelV_given_V.transpose():
                print(*x, file=fw)

            print("# P(delJ | J)", file=fw)
            print(*PdelJ_given_J.transpose().shape, file=fw)
            for x in PdelJ_given_J.transpose():
                print(*x, file=fw)

            print("# P(insVJ)", file=fw)
            print(len(PinsVJ), file=fw)
            print(*PinsVJ, file=fw)

            print("# Markov VJ", file=fw)
            for i in range(4):
                for j in range(4):
                    print(self.generative_model.Rvj[j, i], file=fw, end=" ")
            print(file=fw)

            print("# First nucleotide bias VJ", file=fw)
            print(*first_nt_bias_insVJ, file=fw)

            print("# Error rate", file=fw)
            print(error_rate, file=fw)

            print("# Thymus selection parameter", file=fw)
            print(1., file=fw)

            print("# Conserved J residues", file=fw)
            print("F V W", file=fw)

        return True
