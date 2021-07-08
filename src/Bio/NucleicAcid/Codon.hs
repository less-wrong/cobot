{-# LANGUAGE ViewPatterns #-}

module Bio.NucleicAcid.Codon where

import           Bio.NucleicAcid.Nucleotide (DNA (..), RNA (..), toRNA)
import           Bio.Protein.AminoAcid      (AA (..))
import           Control.Lens               (each, (%~), (&))
import           Control.Monad              ((>=>))
import           Data.Foldable              (Foldable (..))

type Codon a = (a, a, a)

type DNACodon = Codon DNA
type RNACodon = Codon RNA

class Eq a => CodonLike a where
    startCodon :: Codon a
    amberCodon :: Codon a
    ochreCodon :: Codon a
    opalCodon  :: Codon a

    translateCodon :: Codon a -> Maybe AA

    isStartCodon :: Codon a -> Bool
    isStartCodon = (==) startCodon

    isStopCodon :: Codon a -> Bool
    isStopCodon = flip elem [amberCodon, ochreCodon, opalCodon]


instance CodonLike DNA where
    startCodon = (DA, DT, DG)
    amberCodon = (DT, DA, DG)
    ochreCodon = (DT, DA, DA)
    opalCodon  = (DT, DG, DA)

    translateCodon c = c & each %~ toRNA & translateCodon

instance CodonLike RNA where
    startCodon = (RA, RU, RG)
    amberCodon = (RU, RA, RG)
    ochreCodon = (RU, RA, RA)
    opalCodon  = (RU, RG, RA)

    translateCodon = translateRNACodon


splitToTriples :: Foldable t => t a -> Maybe [(a, a, a)]
splitToTriples (toList -> l) | length l `mod` 3 /= 0 = Nothing
                             | otherwise = Just $ helper l
  where
    helper :: [a] -> [(a, a, a)]
    helper []       = []
    helper (x:y:z:xs) = (x, y, z) : helper xs
    helper _        = error "Bio.NucleicAcid.Codon.splitToTriples: impossible"


-- | Translate something with 'CodonLike' to list of aminoacids.
-- Example:
--
-- @
-- ghci> x = "AGAGAGAGAGGAGCGCGCTTAGCTCGATCTTCTATCGAA" :: [DNA]
-- ghci> fmap symbol <$> translateToAA x
-- Just "RERGARLARSSIE"
-- @
--
translateToAA :: (Foldable t, CodonLike a) => t a -> Maybe [AA]
translateToAA = splitToTriples >=> traverse translateCodon



translateRNACodon :: RNACodon -> Maybe AA
translateRNACodon (RU, RU, RU) = Just PHE
translateRNACodon (RU, RU, RC) = Just PHE
translateRNACodon (RU, RU, RA) = Just LEU
translateRNACodon (RU, RU, RG) = Just LEU
translateRNACodon (RC, RU,  _) = Just LEU
translateRNACodon (RA, RU, RG) = Just MET
translateRNACodon (RA, RU,  _) = Just ILE
translateRNACodon (RG, RU,  _) = Just VAL
translateRNACodon (RU, RC,  _) = Just SER
translateRNACodon (RC, RC,  _) = Just PRO
translateRNACodon (RA, RC,  _) = Just THR
translateRNACodon (RG, RC,  _) = Just ALA
translateRNACodon (RU, RA, RU) = Just TYR
translateRNACodon (RU, RA, RC) = Just TYR
translateRNACodon (RU, RA, RA) = Nothing
translateRNACodon (RU, RA, RG) = Nothing
translateRNACodon (RC, RA, RU) = Just HIS
translateRNACodon (RC, RA, RC) = Just HIS
translateRNACodon (RC, RA, RA) = Just GLN
translateRNACodon (RC, RA, RG) = Just GLN
translateRNACodon (RA, RA, RU) = Just ASN
translateRNACodon (RA, RA, RC) = Just ASN
translateRNACodon (RA, RA, RA) = Just LYS
translateRNACodon (RA, RA, RG) = Just LYS
translateRNACodon (RG, RA, RU) = Just ASP
translateRNACodon (RG, RA, RC) = Just ASP
translateRNACodon (RG, RA, RA) = Just GLU
translateRNACodon (RG, RA, RG) = Just GLU
translateRNACodon (RU, RG, RU) = Just CYS
translateRNACodon (RU, RG, RC) = Just CYS
translateRNACodon (RU, RG, RA) = Nothing
translateRNACodon (RU, RG, RG) = Just TRP
translateRNACodon (RC, RG,  _) = Just ARG
translateRNACodon (RA, RG, RU) = Just SER
translateRNACodon (RA, RG, RC) = Just SER
translateRNACodon (RA, RG, RA) = Just ARG
translateRNACodon (RA, RG, RG) = Just ARG
translateRNACodon (RG, RG,  _) = Just GLY
