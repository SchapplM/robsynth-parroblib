% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:50
% EndTime: 2020-03-09 21:14:51
% DurationCPUTime: 0.77s
% Computational Cost: add. (8332->105), mult. (4482->210), div. (1380->5), fcn. (5532->24), ass. (0->126)
t397 = sin(qJ(3,3));
t400 = cos(qJ(3,3));
t391 = pkin(7) + qJ(2,3);
t382 = qJ(3,3) + t391;
t370 = sin(t382);
t373 = cos(t382);
t376 = sin(t391);
t379 = cos(t391);
t458 = 0.1e1 / (t370 * t379 - t373 * t376);
t459 = (mrSges(3,1) * t397 + mrSges(3,2) * t400) * t458;
t392 = pkin(7) + qJ(2,2);
t383 = qJ(3,2) + t392;
t371 = sin(t383);
t374 = cos(t383);
t377 = sin(t392);
t380 = cos(t392);
t457 = 0.1e1 / (t371 * t380 - t374 * t377);
t393 = pkin(7) + qJ(2,1);
t384 = qJ(3,1) + t393;
t372 = sin(t384);
t375 = cos(t384);
t378 = sin(t393);
t381 = cos(t393);
t456 = 0.1e1 / (t372 * t381 - t375 * t378);
t404 = xDP(1);
t408 = 0.1e1 / pkin(2);
t433 = t404 * t408;
t403 = xDP(2);
t434 = t403 * t408;
t395 = legFrame(2,3);
t386 = sin(t395);
t389 = cos(t395);
t361 = -t386 * t371 + t389 * t374;
t438 = t457 * t361;
t360 = t389 * t371 + t386 * t374;
t439 = t457 * t360;
t335 = t433 * t438 + t434 * t439;
t406 = 0.1e1 / pkin(3);
t432 = t406 * t408;
t426 = t457 * t432;
t347 = -pkin(2) * (t377 * t386 - t380 * t389) + t361 * pkin(3);
t443 = t347 * t404;
t344 = pkin(2) * (t377 * t389 + t380 * t386) + t360 * pkin(3);
t446 = t344 * t403;
t410 = (t443 + t446) * t426;
t326 = -t410 + t335;
t455 = pkin(3) * t326;
t396 = legFrame(1,3);
t387 = sin(t396);
t390 = cos(t396);
t363 = -t387 * t372 + t390 * t375;
t436 = t456 * t363;
t362 = t390 * t372 + t387 * t375;
t437 = t456 * t362;
t336 = t433 * t436 + t434 * t437;
t425 = t456 * t432;
t348 = -pkin(2) * (t378 * t387 - t381 * t390) + t363 * pkin(3);
t442 = t348 * t404;
t345 = pkin(2) * (t378 * t390 + t381 * t387) + t362 * pkin(3);
t445 = t345 * t403;
t409 = (t442 + t445) * t425;
t327 = -t409 + t336;
t454 = pkin(3) * t327;
t394 = legFrame(3,3);
t385 = sin(t394);
t388 = cos(t394);
t359 = -t385 * t370 + t388 * t373;
t440 = t458 * t359;
t358 = t388 * t370 + t385 * t373;
t441 = t458 * t358;
t334 = t433 * t440 + t434 * t441;
t427 = t458 * t432;
t346 = -pkin(2) * (t376 * t385 - t379 * t388) + t359 * pkin(3);
t444 = t346 * t404;
t343 = pkin(2) * (t376 * t388 + t379 * t385) + t358 * pkin(3);
t447 = t343 * t403;
t322 = -(t444 / 0.2e1 + t447 / 0.2e1) * t427 + t334;
t411 = (t444 + t447) * t427;
t325 = -t411 + t334;
t405 = pkin(3) ^ 2;
t407 = pkin(2) ^ 2;
t420 = t370 * t376 + t373 * t379;
t414 = t420 * pkin(2);
t431 = 0.2e1 * pkin(2) * pkin(3);
t435 = t458 * t408;
t450 = t325 * t411;
t316 = ((-t420 * t322 * t431 - t325 * t405 - t334 * t407) * t406 * t334 + (pkin(3) + t414) * t450) * t435;
t319 = (-pkin(3) * t450 + (pkin(3) * t325 + t334 * t414) * t334) * t435;
t417 = (-mrSges(3,1) * t400 + mrSges(3,2) * t397) * pkin(2);
t364 = -Ifges(3,3) + t417;
t453 = (Ifges(3,3) * t316 - t319 * t364) * t458;
t323 = -(t443 / 0.2e1 + t446 / 0.2e1) * t426 + t335;
t419 = t371 * t377 + t374 * t380;
t413 = t419 * pkin(2);
t449 = t410 * t457;
t317 = (-(-t419 * t323 * t431 - t326 * t405 - t335 * t407) * t406 * t457 * t335 - t326 * (pkin(3) + t413) * t449) * t408;
t321 = (-t410 * t455 + (t335 * t413 + t455) * t335) * t408 * t457;
t398 = sin(qJ(3,2));
t401 = cos(qJ(3,2));
t416 = (-mrSges(3,1) * t401 + mrSges(3,2) * t398) * pkin(2);
t365 = -Ifges(3,3) + t416;
t452 = (-Ifges(3,3) * t317 - t321 * t365) * t457;
t324 = -(t442 / 0.2e1 + t445 / 0.2e1) * t425 + t336;
t418 = t372 * t378 + t375 * t381;
t412 = t418 * pkin(2);
t448 = t409 * t456;
t318 = (-(-t418 * t324 * t431 - t327 * t405 - t336 * t407) * t406 * t456 * t336 - t327 * (pkin(3) + t412) * t448) * t408;
t320 = (-t409 * t454 + (t336 * t412 + t454) * t336) * t408 * t456;
t399 = sin(qJ(3,1));
t402 = cos(qJ(3,1));
t415 = (-mrSges(3,1) * t402 + mrSges(3,2) * t399) * pkin(2);
t366 = -Ifges(3,3) + t415;
t451 = (-Ifges(3,3) * t318 - t320 * t366) * t456;
t430 = t334 ^ 2 * t459;
t368 = mrSges(3,1) * t398 + mrSges(3,2) * t401;
t429 = t335 ^ 2 * t457 * t368;
t369 = mrSges(3,1) * t399 + mrSges(3,2) * t402;
t428 = t336 ^ 2 * t456 * t369;
t424 = -m(3) * t407 - Ifges(2,3) - Ifges(3,3);
t423 = 0.2e1 * t322 * t411 * t459;
t422 = 0.2e1 * t368 * t323 * t449;
t421 = 0.2e1 * t369 * t324 * t448;
t312 = -(0.2e1 * t415 + t424) * t320 + t366 * t318;
t311 = -(0.2e1 * t416 + t424) * t321 + t365 * t317;
t310 = -(0.2e1 * t417 + t424) * t319 - t364 * t316;
t1 = [t359 * t423 + t361 * t422 + t363 * t421 + (t310 * t440 + t311 * t438 + t312 * t436) * t408 + (-t346 * t430 - t347 * t429 - t348 * t428 + (-t346 * t453 - t347 * t452 - t348 * t451) * t408) * t406; t358 * t423 + t360 * t422 + t362 * t421 + (t310 * t441 + t311 * t439 + t312 * t437) * t408 + (-t343 * t430 - t344 * t429 - t345 * t428 + (-t343 * t453 - t344 * t452 - t345 * t451) * t408) * t406; 0;];
taucX  = t1;
