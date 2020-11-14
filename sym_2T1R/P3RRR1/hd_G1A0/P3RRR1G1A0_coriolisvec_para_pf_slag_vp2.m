% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRR1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:20
% EndTime: 2019-05-03 15:38:22
% DurationCPUTime: 2.08s
% Computational Cost: add. (9424->189), mult. (10634->388), div. (2436->5), fcn. (11542->26), ass. (0->212)
t533 = 0.1e1 / pkin(1);
t501 = qJ(1,2) + qJ(2,2);
t486 = sin(t501);
t489 = cos(t501);
t509 = sin(qJ(1,2));
t515 = cos(qJ(1,2));
t606 = 0.1e1 / (t486 * t515 - t489 * t509);
t585 = t606 * t533;
t500 = qJ(1,3) + qJ(2,3);
t485 = sin(t500);
t488 = cos(t500);
t507 = sin(qJ(1,3));
t513 = cos(qJ(1,3));
t607 = 0.1e1 / (t485 * t513 - t488 * t507);
t588 = t607 * t533;
t608 = 2 * pkin(2);
t502 = qJ(1,1) + qJ(2,1);
t487 = sin(t502);
t490 = cos(t502);
t511 = sin(qJ(1,1));
t517 = cos(qJ(1,1));
t534 = (t487 * t511 + t490 * t517) * pkin(1);
t535 = (t486 * t509 + t489 * t515) * pkin(1);
t536 = (t485 * t507 + t488 * t513) * pkin(1);
t605 = 0.1e1 / (t487 * t517 - t490 * t511);
t518 = xDP(3);
t499 = t518 ^ 2;
t604 = mrSges(2,1) * pkin(1);
t603 = mrSges(2,2) * pkin(1);
t521 = xP(3);
t497 = sin(t521);
t498 = cos(t521);
t524 = koppelP(3,2);
t527 = koppelP(3,1);
t473 = -t497 * t524 + t498 * t527;
t519 = xDP(2);
t446 = t473 * t518 + t519;
t503 = legFrame(3,3);
t491 = sin(t503);
t494 = cos(t503);
t452 = t485 * t494 + t488 * t491;
t440 = pkin(1) * (t491 * t513 + t494 * t507) + t452 * pkin(2);
t531 = 1 / pkin(2);
t578 = t531 * t533;
t560 = t607 * t578;
t549 = t440 * t560;
t428 = t446 * t549;
t470 = t497 * t527 + t498 * t524;
t520 = xDP(1);
t449 = -t470 * t518 + t520;
t453 = -t485 * t491 + t488 * t494;
t443 = pkin(1) * (-t491 * t507 + t494 * t513) + t453 * pkin(2);
t546 = t443 * t560;
t431 = t449 * t546;
t401 = -t431 - t428;
t566 = t453 * t588;
t567 = t452 * t588;
t416 = t446 * t567 + t449 * t566;
t395 = t416 + t401;
t602 = t395 * t401;
t525 = koppelP(2,2);
t528 = koppelP(2,1);
t474 = -t497 * t525 + t498 * t528;
t447 = t474 * t518 + t519;
t504 = legFrame(2,3);
t492 = sin(t504);
t495 = cos(t504);
t454 = t486 * t495 + t489 * t492;
t441 = pkin(1) * (t492 * t515 + t495 * t509) + t454 * pkin(2);
t558 = t606 * t578;
t548 = t441 * t558;
t429 = t447 * t548;
t471 = t497 * t528 + t498 * t525;
t450 = -t471 * t518 + t520;
t455 = -t486 * t492 + t489 * t495;
t444 = pkin(1) * (-t492 * t509 + t495 * t515) + t455 * pkin(2);
t545 = t444 * t558;
t432 = t450 * t545;
t402 = -t432 - t429;
t564 = t455 * t585;
t565 = t454 * t585;
t417 = t447 * t565 + t450 * t564;
t396 = t417 + t402;
t601 = t396 * t402;
t526 = koppelP(1,2);
t529 = koppelP(1,1);
t475 = -t497 * t526 + t498 * t529;
t448 = t475 * t518 + t519;
t505 = legFrame(1,3);
t493 = sin(t505);
t496 = cos(t505);
t456 = t487 * t496 + t490 * t493;
t442 = pkin(1) * (t493 * t517 + t496 * t511) + t456 * pkin(2);
t556 = t605 * t578;
t547 = t442 * t556;
t430 = t448 * t547;
t472 = t497 * t529 + t498 * t526;
t451 = -t472 * t518 + t520;
t457 = -t487 * t493 + t490 * t496;
t445 = pkin(1) * (-t493 * t511 + t496 * t517) + t457 * pkin(2);
t544 = t445 * t556;
t433 = t451 * t544;
t403 = -t433 - t430;
t582 = t605 * t533;
t562 = t457 * t582;
t563 = t456 * t582;
t418 = t448 * t563 + t451 * t562;
t397 = t418 + t403;
t600 = t397 * t403;
t506 = sin(qJ(2,3));
t512 = cos(qJ(2,3));
t479 = mrSges(2,1) * t506 + mrSges(2,2) * t512;
t599 = t416 ^ 2 * t479;
t508 = sin(qJ(2,2));
t514 = cos(qJ(2,2));
t480 = mrSges(2,1) * t508 + mrSges(2,2) * t514;
t598 = t417 ^ 2 * t480;
t510 = sin(qJ(2,1));
t516 = cos(qJ(2,1));
t481 = mrSges(2,1) * t510 + mrSges(2,2) * t516;
t597 = t418 ^ 2 * t481;
t596 = t452 * t607;
t595 = t453 * t607;
t594 = t454 * t606;
t593 = t455 * t606;
t592 = t456 * t605;
t591 = t457 * t605;
t482 = t506 * t603;
t532 = pkin(1) ^ 2;
t543 = m(2) * t532 + Ifges(1,3) + Ifges(2,3);
t576 = t512 * t604;
t467 = -0.2e1 * t482 + t543 + 0.2e1 * t576;
t590 = t607 * t467;
t476 = Ifges(2,3) - t482 + t576;
t589 = t607 * t476;
t483 = t508 * t603;
t575 = t514 * t604;
t468 = -0.2e1 * t483 + t543 + 0.2e1 * t575;
t587 = t606 * t468;
t477 = Ifges(2,3) - t483 + t575;
t586 = t606 * t477;
t484 = t510 * t603;
t574 = t516 * t604;
t469 = -0.2e1 * t484 + t543 + 0.2e1 * t574;
t584 = t605 * t469;
t478 = Ifges(2,3) - t484 + t574;
t583 = t605 * t478;
t581 = t607 * t531;
t580 = t606 * t531;
t579 = t605 * t531;
t577 = t533 * t499;
t573 = t440 * t581;
t572 = t441 * t580;
t571 = t442 * t579;
t570 = t443 * t581;
t569 = t444 * t580;
t568 = t445 * t579;
t561 = t476 * t581;
t559 = t477 * t580;
t557 = t478 * t579;
t392 = -t431 / 0.2e1 - t428 / 0.2e1 + t416;
t555 = -0.2e1 * t392 * t401 * t479;
t393 = -t432 / 0.2e1 - t429 / 0.2e1 + t417;
t554 = -0.2e1 * t393 * t402 * t480;
t394 = -t433 / 0.2e1 - t430 / 0.2e1 + t418;
t553 = -0.2e1 * t394 * t403 * t481;
t552 = t581 * t599;
t551 = t580 * t598;
t550 = t579 * t597;
t542 = t607 * t555;
t541 = t606 * t554;
t540 = t605 * t553;
t530 = pkin(2) ^ 2;
t523 = mrSges(3,1);
t522 = mrSges(3,2);
t427 = (t456 * t475 - t457 * t472) * t582;
t426 = (t454 * t474 - t455 * t471) * t585;
t425 = (t452 * t473 - t453 * t470) * t588;
t424 = (-Ifges(2,3) * t568 + t457 * t583) * t533;
t423 = (-Ifges(2,3) * t571 + t456 * t583) * t533;
t422 = (-Ifges(2,3) * t569 + t455 * t586) * t533;
t421 = (-Ifges(2,3) * t572 + t454 * t586) * t533;
t420 = (-Ifges(2,3) * t570 + t453 * t589) * t533;
t419 = (-Ifges(2,3) * t573 + t452 * t589) * t533;
t412 = (-t445 * t557 + t457 * t584) * t533;
t411 = (-t442 * t557 + t456 * t584) * t533;
t410 = (-t444 * t559 + t455 * t587) * t533;
t409 = (-t441 * t559 + t454 * t587) * t533;
t408 = (-t443 * t561 + t453 * t590) * t533;
t407 = (-t440 * t561 + t452 * t590) * t533;
t406 = (t442 * t475 - t445 * t472) * t556;
t405 = (t441 * t474 - t444 * t471) * t558;
t404 = (t440 * t473 - t443 * t470) * t560;
t400 = -Ifges(2,3) * t406 + t427 * t478;
t399 = -Ifges(2,3) * t405 + t426 * t477;
t398 = -Ifges(2,3) * t404 + t425 * t476;
t391 = -t406 * t478 + t427 * t469;
t390 = -t405 * t477 + t426 * t468;
t389 = -t404 * t476 + t425 * t467;
t388 = (-pkin(2) * t600 + (-t397 * pkin(2) - t418 * t534) * t418) * t582;
t387 = ((-pkin(2) * t396 - t417 * t535) * t417 - pkin(2) * t601) * t585;
t386 = ((-pkin(2) * t395 - t416 * t536) * t416 - pkin(2) * t602) * t588;
t385 = ((t394 * t534 * t608 + t397 * t530 + t532 * t418) * t531 * t418 + (pkin(2) + t534) * t600) * t582;
t384 = ((t393 * t535 * t608 + t396 * t530 + t532 * t417) * t531 * t417 + (pkin(2) + t535) * t601) * t585;
t383 = ((t392 * t536 * t608 + t395 * t530 + t532 * t416) * t531 * t416 + (pkin(2) + t536) * t602) * t588;
t382 = -Ifges(2,3) * t385 - t388 * t478;
t381 = -Ifges(2,3) * t384 - t387 * t477;
t380 = -Ifges(2,3) * t383 - t386 * t476;
t379 = -t385 * t478 - t388 * t469;
t378 = -t384 * t477 - t387 * t468;
t377 = -t383 * t476 - t386 * t467;
t1 = [t379 * t562 - t382 * t544 + t457 * t540 - t445 * t550 + t378 * t564 - t381 * t545 + t455 * t541 - t444 * t551 + t377 * t566 - t380 * t546 + t453 * t542 - t443 * t552 - t499 * (-t497 * t522 + t498 * t523) + (-(t412 * t591 - t424 * t568) * t475 - (t412 * t592 - t424 * t571) * t472 - (t410 * t593 - t422 * t569) * t474 - (t410 * t594 - t422 * t572) * t471 - (t408 * t595 - t420 * t570) * t473 - (t408 * t596 - t420 * t573) * t470) * t577; t379 * t563 - t382 * t547 + t456 * t540 - t442 * t550 + t378 * t565 - t381 * t548 + t454 * t541 - t441 * t551 + t377 * t567 - t380 * t549 + t452 * t542 - t440 * t552 - t499 * (t497 * t523 + t498 * t522) + (-(t411 * t591 - t423 * t568) * t475 - (t411 * t592 - t423 * t571) * t472 - (t409 * t593 - t421 * t569) * t474 - (t409 * t594 - t421 * t572) * t471 - (t407 * t595 - t419 * t570) * t473 - (t407 * t596 - t419 * t573) * t470) * t577; t425 * t377 + t426 * t378 + t427 * t379 - t404 * t380 - t405 * t381 - t406 * t382 + (-(t391 * t591 - t400 * t568) * t475 - (t391 * t592 - t400 * t571) * t472 - (t390 * t593 - t399 * t569) * t474 - (t390 * t594 - t399 * t572) * t471 - (t389 * t595 - t398 * t570) * t473 - (t389 * t596 - t398 * t573) * t470) * t577 + (-t404 * t599 - t405 * t598 - t406 * t597 + t425 * t555 + t426 * t554 + t427 * t553) * pkin(1);];
taucX  = t1;
