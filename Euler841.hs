{-# LANGUAGE BangPatterns #-}

-- Haskell solution.
-- Computes alternating shading area for {p/q} with inradius 1 using a
-- stable cosine-difference formula and parallel evaluation across Fibonacci pairs.

import Control.Concurrent (forkIO)
import Control.Concurrent.Chan (Chan, newChan, readChan, writeChan)
import Control.Monad (forM_, replicateM_, unless)
import Data.Array.IO (IOArray, newArray, getElems, writeArray)
import GHC.Conc (getNumProcessors, setNumCapabilities)
import System.Exit (exitFailure)
import System.IO (hPutStrLn, stderr)
import Text.Printf (printf)

format10 :: Double -> String
format10 x = printf "%.10f" x

checkRounded10 :: String -> Double -> String -> IO ()
checkRounded10 name got expected = do
    let s = format10 got
    unless (s == expected) $ do
        hPutStrLn stderr $ "Validation failed: " ++ name ++ " = " ++ s
        exitFailure

kahanSum :: [Double] -> Double
kahanSum = go 0.0 0.0
  where
    go !s !c [] = s
    go !s !c (x:xs) =
        let !y = x - c
            !t = s + y
            !c' = (t - s) - y
        in go t c' xs

areaA :: Int -> Int -> Double
areaA p q
    | q < 1 = error "q must be >= 1"
    | p <= 2 * q = error "need p > 2q"
    | otherwise =
        let !dp = fromIntegral p :: Double
            !delta = pi / dp
            !sd = sin delta
            !cd = cos delta
            !tanDelta = sd / cd
            go :: Int -> Double -> Double -> Double -> Double -> Double -> Double
            go !k !sinK !cosK !sign !acc !c
                | k > q = acc
                | otherwise =
                    let !sinNext = sinK * cd + cosK * sd
                        !cosNext = cosK * cd - sinK * sd
                        !diff = sd / (cosK * cosNext)
                        !term = sign * diff
                        -- Kahan summation to reduce cancellation.
                        !y = term - c
                        !t = acc + y
                        !c' = (t - acc) - y
                    in go (k + 1) sinNext cosNext (-sign) t c'
            !sumDiff =
                if q < 2
                    then 0.0
                    else go 2 sd cd 1.0 0.0 0.0
            !inner = sumDiff - tanDelta
            !sign = if even q then 1.0 else -1.0
        in dp * sign * inner

fibs :: [Integer]
fibs = 1 : 1 : zipWith (+) fibs (drop 1 fibs)

fibPairs :: [(Int, Int)]
fibPairs =
    let fs = take 36 fibs  -- F1..F36
    in [ (fromIntegral (fs !! n), fromIntegral (fs !! (n - 2)))
       | n <- [3 .. 34] ]

validate :: IO ()
validate = do
    checkRounded10 "A(8,3)" (areaA 8 3) "9.9411254970"
    checkRounded10 "A(130021,50008)" (areaA 130021 50008) "10.9210371479"
    hPutStrLn stderr "Validation checkpoints passed."

worker :: Chan (Maybe (Int, Int, Int)) -> IOArray Int Double -> Chan () -> IO ()
worker taskChan outArr doneChan = do
    task <- readChan taskChan
    case task of
        Nothing -> return ()
        Just (idx, p, q) -> do
            let !val = areaA p q
            writeArray outArr idx val
            writeChan doneChan ()
            worker taskChan outArr doneChan

main :: IO ()
main = do
    procs <- getNumProcessors
    let workers = max 1 procs
    setNumCapabilities workers

    validate

    let pairs = fibPairs
        n = length pairs
    outArr <- newArray (0, n - 1) 0.0 :: IO (IOArray Int Double)
    taskChan <- newChan
    doneChan <- newChan

    replicateM_ workers (forkIO (worker taskChan outArr doneChan))

    forM_ (zip [0 ..] pairs) $ \(idx, (p, q)) ->
        writeChan taskChan (Just (idx, p, q))
    replicateM_ workers (writeChan taskChan Nothing)

    replicateM_ n (readChan doneChan)

    vals <- getElems outArr
    let total = kahanSum vals
    printf "%.10f\n" total
