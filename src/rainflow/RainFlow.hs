{- |
Module      :  $Header$
Description :  Simple implementation of RainFlow cycle count algorithm for fatigue analysis (as ASTM E1049 − 85)
Copyright   :  (c) 2018 Fabio 'Monte' Sant'Anna
License     :  <license>

Maintainer  :  Fabio 'Monte' Sant'Anna <fsantanna@zoho.com>
Stability   :  provisional
Portability :  portable

Simple implementation of RainFlow cycle count algorithm for fatigue analysis, as
a direct adaptation of 5.4.4 algorithm ASTM E1049-85 [1].

Also adequate data types and functions are provided to manipulate:
 - stress histories;
 - simple alternate / mean stresses pair cycle counts;
 - bin alternate / mean stresses pair cycle counts.

[1] Standard, A. S. T. M. "E1049-85." Standard Practices for Cycle Counting in Fatigue Analysis, ASTM (2011).
-}

-- {-# LANGUAGE  GeneralizedNewtypeDeriving  #-}

module RainFlow
  ( Stress, stress
  , Stresses, emptyStresses, stressesFromList, filterPeaks, trimExtremities
  , CycleCount, zeroCycle, halfCycle, oneCycle, cycles
  , RFRanges, emptyRFRanges, rainflow
  , RFBins, newRFBinCounts, insRangeBin, rangesToBins)
  where

import Data.Monoid
-- import Data.List (intersperse)
-- import Data.Foldable
import qualified Data.Vector.Unboxed as V
import qualified Data.Map.Strict as M


-- | A stress value type
type Stress = Double

-- | Constructs a stress value from a Double
stress :: Double -> Stress
stress s = s

-- | An array of stresses
newtype Stresses = Sts (V.Vector Stress) deriving (Ord, Eq)

instance Show Stresses where
  show (Sts xs) = show $ V.toList xs

-- | Constructs an empty Stress array
emptyStresses :: Stresses
emptyStresses = Sts V.empty -- :: V.Vector Stress

-- | Constructs an Stress array from a list
stressesFromList :: [Double] -> Stresses
stressesFromList = Sts . V.fromList

-- | Leave only stress peaks and valleys in the stress array
filterPeaks :: Stresses -> Stresses
filterPeaks (Sts xs)
  | V.length xs <=2 = Sts xs
  | otherwise = let
      ss = fst $ V.foldl' fPeaks (V.empty :: V.Vector Stress, V.head xs) (V.tail xs)
      lx = V.last xs in
      Sts (V.reverse $ if V.head ss /= lx then lx `V.cons` ss else ss)
    where
      fPeaks :: (V.Vector Stress, Stress) -> Stress -> (V.Vector Stress, Stress)
      fPeaks (ys, y) x
        | V.null ys = (V.singleton y, x)
        | (y' > y && x > y) || (y' < y && x < y) = (V.cons y ys, x)
        | otherwise = (ys, x)
          where
            y' = V.head ys

-- | Trim the head and the last stresses. Useful in case of regular cycles with
-- irrelevant data as begin / end of stresses list.
trimExtremities :: Stresses -> Stresses
trimExtremities (Sts xs)
  | V.length xs <= 2 = emptyStresses
  | otherwise = Sts . V.tail . V.init $ xs

{-- Full (no adittional Half Cycle)
    or Half cycle (only apply to more advance cycle count techniques) --}
data HalfCycle = NoHC | HC deriving (Show, Ord, Eq)

-- | Cycles counted. With a resolution of 0.5 cycles.
data CycleCount = CyC
    { cCount :: Word
    , hCycle :: HalfCycle} deriving (Ord, Eq)

instance Show CycleCount where
  show cc = show (fromIntegral (cCount cc) + (if hCycle cc == HC then 0.5 else 0) :: Float)

-- Monoid CycleCount
instance Monoid CycleCount where
  mempty = zeroCycle
  mappend = sumCycles

-- | Cycle constructor: 0.0 cycle
zeroCycle :: CycleCount
zeroCycle = CyC 0 NoHC

-- | Cycle constructor: 0.5 cycle
halfCycle :: CycleCount
halfCycle = CyC 0 HC

-- | Cycle constructor: 1.0 cycle
oneCycle :: CycleCount
oneCycle = CyC 1 NoHC

-- | Cycle constructor: n cycles
cycles :: Word -> CycleCount
cycles n = CyC n NoHC

sumCycles :: CycleCount -> CycleCount -> CycleCount
sumCycles (CyC x NoHC) (CyC y NoHC)  = CyC (x + y) NoHC
sumCycles (CyC x HC)   (CyC y NoHC)  = CyC (x + y) HC
sumCycles (CyC x NoHC) (CyC y HC)    = CyC (x + y) HC
sumCycles (CyC x HC)   (CyC y HC)    = CyC (x + y + 1) NoHC


-- | Represents a count of alternate / mean stresses pair result of a rainflow (or other
--stress count thechnique) witout bin aggregation
newtype RFRanges = RFRs (M.Map (Stress, Stress) CycleCount) deriving (Ord, Eq)
-- type RFRanges = M.Map (Stress, Stress) CycleCount

instance Show RFRanges where
  show (RFRs m) = M.foldlWithKey' (\ ls (aSt, mSt) cc ->
    concat [ls, " ", show aSt, " ", show mSt, " ", show cc, "\n"]) "" m
    -- concat $ intersperse " " [show ls, show aSt, show mSt] ++ ["\n"]) "" m
  -- ls ++ show k ++ show cc ++ "\n") "" m

-- | Constructs an empty block of rainflow counts
emptyRFRanges :: RFRanges
emptyRFRanges = RFRs M.empty

-- | Inserts a range in a block of rainflow counts. If the range whith the
-- specified alternate / mean stress exist, the count number provided is added
-- to existing in the block.
insRange :: Stress -- ^ Alternate stress
         -> Stress -- ^ Mean stress
         -> CycleCount  -- ^ Cycle count to sum
         -> RFRanges -- ^ Rainflow count data
         -> RFRanges-- ^ new rainflow count data
insRange aSt mSt cc (RFRs m) = RFRs (M.insertWith (<>) (aSt, mSt) cc m)

-- | Performs the basic rainflow algorithm in an array of stresses. The stresses
--should be already filtered (see filterPeaks) or strange things will happen.
-- The output will be NOT aggregated in bins.
rainflow :: Stresses -- ^ Stress 'list'
         -> RFRanges -- ^ Counted cycles in stress blocks
rainflow (Sts xs)
  | V.length xs <= 2 = halfCycles xs emptyRFRanges
  | otherwise = rainflow' emptyRFRanges (V.tail xs) (V.head xs) (V.empty :: V.Vector Stress)
    where
      rainflow' :: RFRanges -> V.Vector Stress -> Stress -> V.Vector Stress -> RFRanges
      rainflow' acc zs p ys
        | V.length zs == 0 = halfCycles (V.cons p ys) acc
        | V.length ys == 0 = rainflow' acc (V.tail zs) (V.head zs) (V.cons p ys)
        | otherwise         =
          let x = V.head zs; y = V.head ys; rX = abs (x-p); rY = abs (y-p) in
          if rX < rY then rainflow' acc (V.tail zs) x (V.cons p ys) else
            if V.length ys == 1 then
              rainflow' (insRange rY ((y+p)/2) halfCycle acc)
                      (V.tail zs) x (V.cons p $ V.tail ys)  -- step (5)
            else rainflow' (insRange rY ((y+p)/2) oneCycle acc)
                        (V.tail zs) x (V.tail ys)  -- step (4)

      halfCycles :: V.Vector Stress -> RFRanges -> RFRanges
      halfCycles ys rfR
        | V.length ys <= 1 = rfR
        | otherwise = fst $ V.foldl' (\(cc, x) x' ->
              (insRange (abs (x-x')) ((x+x')/2) halfCycle cc, x')) (rfR, y) ys'
        where
          y = V.head ys; ys' = V.tail ys


-- | Represents a count of alternate / mean stresses pair result of a rainflow (or other
--stress count thechnique) in bin aggregation with a defined stress
data RFBins = RFBs { aBinSize :: Stress, mBinSize :: Stress
                         , ranges :: M.Map (Int, Int) CycleCount }
    deriving (Show, Eq)

-- | Constructs an empty block of rainflow counts, with provided stress as bins
newRFBinCounts :: Stress -- ^ Alternate stress
               -> Stress -- ^ Mean stress
               -> RFBins -- ^ New empty block of rainflow counts
newRFBinCounts binSizeA binSizeM = RFBs binSizeA binSizeM M.empty


binFromStress :: Stress -> Stress -> Int
binFromStress binSize st = truncate $ st / binSize

stressFromBin :: Stress -> Int -> Stress
stressFromBin binSize bin = fromIntegral bin * binSize

stressesFromBin :: Stress -> Int -> (Stress, Stress)
stressesFromBin binSize bin = (minStBin, minStBin + binSize)
  where minStBin = stressFromBin binSize bin

-- | Inserts a range in a block of rainflow bin counts. If the range whith the
-- specified alternate / mean stress exist, the count number provided is added
-- to existing in the block.
insRangeBin :: Stress -- ^ Alternate stress
            -> Stress -- ^ Mean stress
            -> CycleCount -- ^ Cycle count to sum
            -> RFBins -- ^ Rainflow bin count data
            -> RFBins -- ^ New rainflow bin count data
insRangeBin aSt mSt cc rfb =
   RFBs binSzA binSzM (M.insertWith (<>) (binA, binM) cc (ranges rfb))
    where
      binSzA = aBinSize rfb
      binSzM = mBinSize rfb
      binA = binFromStress binSzA aSt
      binM = binFromStress binSzM mSt

-- | Adjusts a simple RFRanges block ou count data in a bin configuration. If
-- there is data in the RFBins variable, the new cycles are added.
rangesToBins :: RFRanges -- ^ Rainflow count data
             -> RFBins -- ^ Rainflow bin count data
             -> RFBins -- ^ New rainflow bin count data
rangesToBins (RFRs rfR) rfB = M.foldrWithKey (\(aSt, mSt) cc rfb ->
    insRangeBin aSt mSt cc rfb) rfB rfR
