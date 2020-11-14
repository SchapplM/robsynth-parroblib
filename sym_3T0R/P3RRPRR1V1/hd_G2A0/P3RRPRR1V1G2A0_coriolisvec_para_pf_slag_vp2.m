% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
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
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:33:54
% EndTime: 2020-08-06 19:33:55
% DurationCPUTime: 1.14s
% Computational Cost: add. (5292->171), mult. (5799->333), div. (2331->7), fcn. (5733->18), ass. (0->145)
t528 = (pkin(1) * m(3));
t527 = -2 * mrSges(3,3);
t453 = legFrame(3,2);
t432 = sin(t453);
t435 = cos(t453);
t462 = cos(qJ(2,3));
t443 = 0.1e1 / t462;
t471 = pkin(2) + pkin(1);
t449 = 0.1e1 / t471;
t469 = xDP(2);
t470 = xDP(1);
t402 = (t432 * t470 + t435 * t469) * t449 * t443;
t399 = t402 ^ 2;
t454 = legFrame(2,2);
t433 = sin(t454);
t436 = cos(t454);
t464 = cos(qJ(2,2));
t445 = 0.1e1 / t464;
t403 = (t433 * t470 + t436 * t469) * t449 * t445;
t400 = t403 ^ 2;
t455 = legFrame(1,2);
t434 = sin(t455);
t437 = cos(t455);
t466 = cos(qJ(2,1));
t447 = 0.1e1 / t466;
t404 = (t434 * t470 + t437 * t469) * t449 * t447;
t401 = t404 ^ 2;
t448 = t471 ^ 2;
t456 = sin(qJ(2,3));
t457 = sin(qJ(1,3));
t507 = t457 * t462;
t411 = -t432 * t507 + t456 * t435;
t414 = t456 * t432 + t435 * t507;
t450 = pkin(3) + qJ(3,3);
t439 = 0.1e1 / t450;
t463 = cos(qJ(1,3));
t468 = xDP(3);
t384 = (t463 * t468 + (t411 * t469 + t414 * t470) * t443) * t439;
t526 = 0.2e1 * t384;
t458 = sin(qJ(2,2));
t459 = sin(qJ(1,2));
t504 = t459 * t464;
t412 = -t433 * t504 + t458 * t436;
t415 = t458 * t433 + t436 * t504;
t451 = pkin(3) + qJ(3,2);
t440 = 0.1e1 / t451;
t465 = cos(qJ(1,2));
t385 = (t465 * t468 + (t412 * t469 + t415 * t470) * t445) * t440;
t525 = 0.2e1 * t385;
t460 = sin(qJ(2,1));
t461 = sin(qJ(1,1));
t501 = t461 * t466;
t413 = -t434 * t501 + t460 * t437;
t416 = t460 * t434 + t437 * t501;
t452 = pkin(3) + qJ(3,1);
t441 = 0.1e1 / t452;
t467 = cos(qJ(1,1));
t386 = (t467 * t468 + (t413 * t469 + t416 * t470) * t447) * t441;
t524 = 0.2e1 * t386;
t482 = (-2 * mrSges(3,1) - t528) * pkin(1);
t518 = (-Ifges(2,1) - Ifges(3,1));
t423 = Ifges(2,2) + Ifges(3,2) - t482 + t518;
t523 = 2 * t423;
t428 = mrSges(3,2) * pkin(1) - Ifges(2,4) - Ifges(3,4);
t522 = 2 * t428;
t521 = (m(3) * qJ(3,1));
t520 = (m(3) * qJ(3,2));
t519 = (m(3) * qJ(3,3));
t517 = (-Ifges(2,5) - Ifges(3,5));
t516 = (-Ifges(2,6) - Ifges(3,6));
t515 = t399 * t456;
t514 = t400 * t458;
t513 = t401 * t460;
t442 = t462 ^ 2;
t512 = t428 * t442;
t444 = t464 ^ 2;
t511 = t428 * t444;
t446 = t466 ^ 2;
t510 = t428 * t446;
t509 = t456 * t462;
t508 = t456 * t471;
t506 = t458 * t464;
t505 = t458 * t471;
t503 = t460 * t466;
t502 = t460 * t471;
t500 = t471 * t462;
t499 = t471 * t464;
t498 = t471 * t466;
t475 = -t450 * t463 + t457 * t500;
t393 = -t475 * t432 + t435 * t508;
t396 = t432 * t508 + t475 * t435;
t417 = t457 * t450 + t463 * t500;
t381 = (t393 * t469 + t396 * t470 + t417 * t468) * t439;
t494 = 0.2e1 * t471;
t363 = (-t448 * t399 + ((-t448 * t442 - t450 ^ 2) * t384 + (t402 * t456 * t450 + t381 * t462) * t494) * t384) * t439;
t372 = (-t471 * t399 * t443 + (-t384 * t500 + 0.2e1 * t381) * t384) * t439;
t425 = mrSges(3,2) * qJ(3,3) + t516;
t438 = mrSges(3,1) + t528;
t483 = -pkin(1) * mrSges(3,3) - t517;
t405 = -(-t438 * qJ(3,3) + t483) * t456 + t462 * t425;
t420 = -t456 * mrSges(3,2) + t438 * t462;
t429 = mrSges(3,3) + t519;
t492 = t443 * t515;
t493 = -Ifges(1,3) + t518;
t497 = (t509 * t522 - t423 * t442 + ((t527 - t519) * qJ(3,3)) + t493) * t372 - t405 * t492 + t420 * t363 - 0.4e1 * t384 * t402 * t512 - (t384 * t456 * t523 + (qJ(3,3) * mrSges(3,1) + t429 * pkin(1) + t517) * t402) * t402 * t462 + t425 * t515 + (t381 * t429 + t428 * t402) * t526;
t474 = -t451 * t465 + t459 * t499;
t394 = -t474 * t433 + t436 * t505;
t397 = t433 * t505 + t474 * t436;
t418 = t459 * t451 + t465 * t499;
t382 = (t394 * t469 + t397 * t470 + t418 * t468) * t440;
t364 = (-t448 * t400 + ((-t448 * t444 - t451 ^ 2) * t385 + (t403 * t458 * t451 + t382 * t464) * t494) * t385) * t440;
t373 = (-t471 * t400 * t445 + (-t385 * t499 + 0.2e1 * t382) * t385) * t440;
t426 = mrSges(3,2) * qJ(3,2) + t516;
t406 = -(-t438 * qJ(3,2) + t483) * t458 + t464 * t426;
t421 = -t458 * mrSges(3,2) + t438 * t464;
t430 = mrSges(3,3) + t520;
t491 = t445 * t514;
t496 = (t506 * t522 - t423 * t444 + ((t527 - t520) * qJ(3,2)) + t493) * t373 - t406 * t491 + t421 * t364 - 0.4e1 * t385 * t403 * t511 - (t385 * t458 * t523 + (qJ(3,2) * mrSges(3,1) + t430 * pkin(1) + t517) * t403) * t403 * t464 + t426 * t514 + (t382 * t430 + t428 * t403) * t525;
t473 = -t452 * t467 + t461 * t498;
t395 = -t473 * t434 + t437 * t502;
t398 = t434 * t502 + t473 * t437;
t419 = t461 * t452 + t467 * t498;
t383 = (t395 * t469 + t398 * t470 + t419 * t468) * t441;
t365 = (-t448 * t401 + ((-t448 * t446 - t452 ^ 2) * t386 + (t404 * t460 * t452 + t383 * t466) * t494) * t386) * t441;
t374 = (-t471 * t401 * t447 + (-t386 * t498 + 0.2e1 * t383) * t386) * t441;
t427 = mrSges(3,2) * qJ(3,1) + t516;
t407 = -(-t438 * qJ(3,1) + t483) * t460 + t466 * t427;
t422 = -t460 * mrSges(3,2) + t438 * t466;
t431 = mrSges(3,3) + t521;
t490 = t447 * t513;
t495 = (t503 * t522 - t423 * t446 + ((t527 - t521) * qJ(3,1)) + t493) * t374 - t407 * t490 + t422 * t365 - 0.4e1 * t386 * t404 * t510 - (t386 * t460 * t523 + (qJ(3,1) * mrSges(3,1) + t431 * pkin(1) + t517) * t404) * t404 * t466 + t427 * t513 + (t383 * t431 + t428 * t404) * t524;
t489 = t497 * t443;
t488 = t496 * t445;
t487 = t495 * t447;
t481 = t462 * mrSges(3,2) + t438 * t456;
t486 = -m(3) * t363 + t420 * t372 + (-t384 * t429 / 0.2e1 + t481 * t402) * t526;
t480 = t464 * mrSges(3,2) + t438 * t458;
t485 = -m(3) * t364 + t421 * t373 + (-t385 * t430 / 0.2e1 + t480 * t403) * t525;
t479 = t466 * mrSges(3,2) + t438 * t460;
t484 = -m(3) * t365 + t422 * t374 + (-t386 * t431 / 0.2e1 + t479 * t404) * t524;
t424 = -Ifges(2,3) - Ifges(3,3) + t482;
t478 = (t405 * t372 + (-0.2e1 * t481 * t381 + (t423 * t509 - t428 + 0.2e1 * t512) * t384) * t384 - t424 * t492) * t443;
t477 = (t406 * t373 + (-0.2e1 * t480 * t382 + (t423 * t506 - t428 + 0.2e1 * t511) * t385) * t385 - t424 * t491) * t445;
t476 = (t407 * t374 + (-0.2e1 * t479 * t383 + (t423 * t503 - t428 + 0.2e1 * t510) * t386) * t386 - t424 * t490) * t447;
t1 = [(t484 * t398 + t416 * t487) * t441 + (t485 * t397 + t415 * t488) * t440 + (t486 * t396 + t414 * t489) * t439 + (t432 * t478 + t433 * t477 + t434 * t476) * t449; (t484 * t395 + t413 * t487) * t441 + (t485 * t394 + t412 * t488) * t440 + (t486 * t393 + t411 * t489) * t439 + (t435 * t478 + t436 * t477 + t437 * t476) * t449; (t484 * t419 + t495 * t467) * t441 + (t485 * t418 + t496 * t465) * t440 + (t486 * t417 + t497 * t463) * t439;];
taucX  = t1;
