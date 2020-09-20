% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR2G2P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:35
% EndTime: 2020-03-09 21:21:36
% DurationCPUTime: 0.83s
% Computational Cost: add. (2454->116), mult. (4590->225), div. (2610->5), fcn. (5070->21), ass. (0->141)
t385 = legFrame(1,2);
t379 = cos(t385);
t399 = xDP(2);
t445 = t379 * t399;
t376 = sin(t385);
t400 = xDP(1);
t448 = t376 * t400;
t474 = t445 + t448;
t384 = legFrame(2,2);
t378 = cos(t384);
t446 = t378 * t399;
t375 = sin(t384);
t449 = t375 * t400;
t473 = t446 + t449;
t383 = legFrame(3,2);
t377 = cos(t383);
t447 = t377 * t399;
t374 = sin(t383);
t450 = t374 * t400;
t472 = t447 + t450;
t392 = cos(qJ(3,3));
t471 = pkin(1) * t392;
t394 = cos(qJ(3,2));
t470 = pkin(1) * t394;
t396 = cos(qJ(3,1));
t469 = pkin(1) * t396;
t402 = 0.1e1 / pkin(2);
t398 = xDP(3);
t404 = 0.1e1 / pkin(1);
t438 = t398 * t404;
t421 = t402 * t438;
t371 = sin(qJ(2,3) + qJ(3,3));
t387 = sin(qJ(2,3));
t365 = pkin(1) * t387 + pkin(2) * t371;
t386 = sin(qJ(3,3));
t380 = 0.1e1 / t386;
t462 = t365 * t380;
t353 = t421 * t462;
t393 = cos(qJ(2,3));
t441 = t386 * t387;
t356 = (pkin(2) * t392 + pkin(1)) * t393 - pkin(2) * t441;
t465 = t356 * t402;
t407 = t472 * t465;
t359 = t392 * t393 - t441;
t444 = t380 * t404;
t436 = t472 * t359 * t444;
t455 = t371 * t398;
t335 = t353 + (-t407 - t455) * t444 + t436;
t338 = -t407 * t444 + t353;
t468 = t335 * t338;
t372 = sin(qJ(2,2) + qJ(3,2));
t389 = sin(qJ(2,2));
t366 = pkin(1) * t389 + pkin(2) * t372;
t388 = sin(qJ(3,2));
t381 = 0.1e1 / t388;
t461 = t366 * t381;
t354 = t421 * t461;
t395 = cos(qJ(2,2));
t440 = t388 * t389;
t357 = (pkin(2) * t394 + pkin(1)) * t395 - pkin(2) * t440;
t464 = t357 * t402;
t406 = t473 * t464;
t360 = t394 * t395 - t440;
t443 = t381 * t404;
t435 = t473 * t360 * t443;
t453 = t372 * t398;
t336 = t354 + (-t406 - t453) * t443 + t435;
t339 = -t406 * t443 + t354;
t467 = t336 * t339;
t373 = sin(qJ(2,1) + qJ(3,1));
t391 = sin(qJ(2,1));
t367 = pkin(1) * t391 + pkin(2) * t373;
t390 = sin(qJ(3,1));
t382 = 0.1e1 / t390;
t460 = t367 * t382;
t355 = t421 * t460;
t397 = cos(qJ(2,1));
t439 = t390 * t391;
t358 = (pkin(2) * t396 + pkin(1)) * t397 - pkin(2) * t439;
t463 = t358 * t402;
t405 = t474 * t463;
t361 = t396 * t397 - t439;
t442 = t382 * t404;
t434 = t474 * t361 * t442;
t451 = t373 * t398;
t337 = t355 + (-t405 - t451) * t442 + t434;
t340 = -t405 * t442 + t355;
t466 = t337 * t340;
t459 = (mrSges(3,1) * t386 + mrSges(3,2) * t392) * t380;
t458 = (mrSges(3,1) * t388 + mrSges(3,2) * t394) * t381;
t457 = (mrSges(3,1) * t390 + mrSges(3,2) * t396) * t382;
t456 = t371 * t380;
t454 = t372 * t381;
t452 = t373 * t382;
t437 = 0.2e1 * pkin(1) * pkin(2);
t332 = t353 / 0.2e1 + (-t455 + (-t450 / 0.2e1 - t447 / 0.2e1) * t465) * t444 + t436;
t344 = -t438 * t456 + t436;
t401 = pkin(2) ^ 2;
t403 = pkin(1) ^ 2;
t327 = ((t332 * t392 * t437 + t335 * t401 + t344 * t403) * t402 * t344 + (pkin(2) + t471) * t468) * t444;
t329 = ((-pkin(2) * t335 - t344 * t471) * t344 - pkin(2) * t468) * t444;
t410 = (-mrSges(3,1) * t392 + mrSges(3,2) * t386) * pkin(1);
t362 = -Ifges(3,3) + t410;
t414 = -m(3) * t403 - Ifges(2,3) - Ifges(3,3);
t320 = (0.2e1 * t410 + t414) * t329 + t362 * t327;
t433 = t320 * t359 * t380;
t333 = t354 / 0.2e1 + (-t453 + (-t449 / 0.2e1 - t446 / 0.2e1) * t464) * t443 + t435;
t345 = -t438 * t454 + t435;
t328 = ((t333 * t394 * t437 + t336 * t401 + t345 * t403) * t402 * t345 + (pkin(2) + t470) * t467) * t443;
t330 = ((-pkin(2) * t336 - t345 * t470) * t345 - pkin(2) * t467) * t443;
t409 = (-mrSges(3,1) * t394 + mrSges(3,2) * t388) * pkin(1);
t363 = -Ifges(3,3) + t409;
t321 = (0.2e1 * t409 + t414) * t330 + t363 * t328;
t432 = t321 * t360 * t381;
t334 = t355 / 0.2e1 + (-t451 + (-t448 / 0.2e1 - t445 / 0.2e1) * t463) * t442 + t434;
t346 = -t438 * t452 + t434;
t326 = ((t334 * t396 * t437 + t337 * t401 + t346 * t403) * t402 * t346 + (pkin(2) + t469) * t466) * t442;
t331 = ((-pkin(2) * t337 - t346 * t469) * t346 - pkin(2) * t466) * t442;
t408 = (-mrSges(3,1) * t396 + mrSges(3,2) * t390) * pkin(1);
t364 = -Ifges(3,3) + t408;
t322 = (0.2e1 * t408 + t414) * t331 + t364 * t326;
t431 = t322 * t361 * t382;
t323 = -Ifges(3,3) * t327 + t329 * t362;
t430 = t323 * t356 * t380;
t324 = -Ifges(3,3) * t328 + t330 * t363;
t429 = t324 * t357 * t381;
t325 = -Ifges(3,3) * t326 + t331 * t364;
t428 = t325 * t358 * t382;
t427 = t344 ^ 2 * t459;
t426 = t345 ^ 2 * t458;
t425 = t346 ^ 2 * t457;
t420 = t332 * t338 * t459;
t419 = t333 * t339 * t458;
t418 = t334 * t340 * t457;
t417 = t356 * t427;
t416 = t357 * t426;
t415 = t358 * t425;
t413 = -0.2e1 * t359 * t420;
t412 = -0.2e1 * t360 * t419;
t411 = -0.2e1 * t361 * t418;
t1 = [t374 * t413 + t375 * t412 + t376 * t411 + (t374 * t433 + t375 * t432 + t376 * t431) * t404 + (-t374 * t417 - t375 * t416 - t376 * t415 + (-t374 * t430 - t375 * t429 - t376 * t428) * t404) * t402; t377 * t413 + t378 * t412 + t379 * t411 + (t377 * t433 + t378 * t432 + t379 * t431) * t404 + (-t377 * t417 - t378 * t416 - t379 * t415 + (-t377 * t430 - t378 * t429 - t379 * t428) * t404) * t402; 0.2e1 * t371 * t420 + 0.2e1 * t372 * t419 + 0.2e1 * t373 * t418 + (-t320 * t456 - t321 * t454 - t322 * t452) * t404 + (t365 * t427 + t366 * t426 + t367 * t425 + (t323 * t462 + t324 * t461 + t325 * t460) * t404) * t402;];
taucX  = t1;
