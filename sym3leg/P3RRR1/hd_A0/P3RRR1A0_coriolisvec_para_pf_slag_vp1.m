% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRR1A0
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
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
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

function taucX = P3RRR1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:15
% EndTime: 2019-05-03 15:38:17
% DurationCPUTime: 2.20s
% Computational Cost: add. (10090->199), mult. (11590->398), div. (2436->5), fcn. (11938->20), ass. (0->193)
t537 = 0.1e1 / pkin(1);
t509 = qJ(1,2) + qJ(2,2);
t489 = sin(t509);
t492 = cos(t509);
t515 = sin(qJ(1,2));
t518 = cos(qJ(1,2));
t604 = 0.1e1 / (t489 * t518 - t515 * t492);
t588 = t604 * t537;
t508 = qJ(1,3) + qJ(2,3);
t488 = sin(t508);
t491 = cos(t508);
t514 = sin(qJ(1,3));
t517 = cos(qJ(1,3));
t605 = 0.1e1 / (t488 * t517 - t514 * t491);
t589 = t605 * t537;
t523 = xP(3);
t505 = sin(t523);
t506 = cos(t523);
t526 = koppelP(3,2);
t529 = koppelP(3,1);
t477 = -t505 * t526 + t506 * t529;
t520 = xDP(3);
t521 = xDP(2);
t453 = t477 * t520 + t521;
t474 = t505 * t529 + t506 * t526;
t522 = xDP(1);
t456 = -t474 * t520 + t522;
t511 = legFrame(3,3);
t494 = sin(t511);
t497 = cos(t511);
t460 = -t488 * t494 + t497 * t491;
t447 = pkin(1) * (-t514 * t494 + t497 * t517) + t460 * pkin(2);
t535 = 0.1e1 / pkin(2);
t582 = t535 * t537;
t570 = t605 * t582;
t559 = t447 * t570;
t459 = t497 * t488 + t494 * t491;
t444 = pkin(1) * (t494 * t517 + t497 * t514) + t459 * pkin(2);
t562 = t444 * t570;
t405 = -t453 * t562 - t456 * t559;
t527 = koppelP(2,2);
t530 = koppelP(2,1);
t478 = -t505 * t527 + t506 * t530;
t454 = t478 * t520 + t521;
t475 = t505 * t530 + t506 * t527;
t457 = -t475 * t520 + t522;
t512 = legFrame(2,3);
t495 = sin(t512);
t498 = cos(t512);
t462 = -t489 * t495 + t498 * t492;
t448 = pkin(1) * (-t515 * t495 + t498 * t518) + t462 * pkin(2);
t568 = t604 * t582;
t558 = t448 * t568;
t461 = t498 * t489 + t495 * t492;
t445 = pkin(1) * (t495 * t518 + t498 * t515) + t461 * pkin(2);
t561 = t445 * t568;
t406 = -t454 * t561 - t457 * t558;
t528 = koppelP(1,2);
t531 = koppelP(1,1);
t479 = -t505 * t528 + t506 * t531;
t455 = t479 * t520 + t521;
t476 = t505 * t531 + t506 * t528;
t458 = -t476 * t520 + t522;
t510 = qJ(1,1) + qJ(2,1);
t490 = sin(t510);
t493 = cos(t510);
t513 = legFrame(1,3);
t496 = sin(t513);
t499 = cos(t513);
t464 = -t490 * t496 + t499 * t493;
t516 = sin(qJ(1,1));
t519 = cos(qJ(1,1));
t449 = pkin(1) * (-t516 * t496 + t499 * t519) + t464 * pkin(2);
t603 = 0.1e1 / (t490 * t519 - t516 * t493);
t566 = t603 * t582;
t557 = t449 * t566;
t463 = t499 * t490 + t496 * t493;
t446 = pkin(1) * (t496 * t519 + t499 * t516) + t463 * pkin(2);
t560 = t446 * t566;
t407 = -t455 * t560 - t458 * t557;
t420 = (t453 * t459 + t456 * t460) * t589;
t611 = 0.2e1 * t420 + t405;
t421 = (t454 * t461 + t457 * t462) * t588;
t610 = 0.2e1 * t421 + t406;
t587 = t603 * t537;
t422 = (t455 * t463 + t458 * t464) * t587;
t609 = 0.2e1 * t422 + t407;
t538 = (t490 * t516 + t493 * t519) * pkin(1);
t539 = (t489 * t515 + t492 * t518) * pkin(1);
t540 = (t488 * t514 + t491 * t517) * pkin(1);
t398 = t422 + t407;
t534 = pkin(2) ^ 2;
t536 = pkin(1) ^ 2;
t596 = t398 * t407;
t383 = ((pkin(2) * t538 * t609 + t398 * t534 + t536 * t422) * t535 * t422 + (pkin(2) + t538) * t596) * t587;
t386 = (-pkin(2) * t596 + (-pkin(2) * t398 - t422 * t538) * t422) * t587;
t503 = rSges(2,1) ^ 2 + rSges(2,2) ^ 2;
t487 = t536 + t503;
t482 = -rSges(2,1) * t516 + rSges(2,2) * t519;
t485 = rSges(2,1) * t519 + rSges(2,2) * t516;
t541 = (-t482 * t490 + t485 * t493) * pkin(1);
t553 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(1,3) + Icges(2,3);
t440 = (t487 + 0.2e1 * t541) * m(2) + t553;
t443 = Icges(2,3) + (t503 + t541) * m(2);
t377 = -t383 * t443 - t386 * t440;
t600 = m(2) * (-t482 * t493 - t485 * t490);
t554 = t407 * t600 * t609;
t608 = t377 * t587 + t603 * t554;
t397 = t421 + t406;
t597 = t397 * t406;
t382 = ((pkin(2) * t539 * t610 + t397 * t534 + t536 * t421) * t535 * t421 + (pkin(2) + t539) * t597) * t588;
t385 = ((-pkin(2) * t397 - t421 * t539) * t421 - pkin(2) * t597) * t588;
t481 = -rSges(2,1) * t515 + rSges(2,2) * t518;
t484 = rSges(2,1) * t518 + rSges(2,2) * t515;
t542 = (-t481 * t489 + t484 * t492) * pkin(1);
t439 = (t487 + 0.2e1 * t542) * m(2) + t553;
t442 = Icges(2,3) + (t503 + t542) * m(2);
t376 = -t382 * t442 - t385 * t439;
t601 = m(2) * (-t481 * t492 - t484 * t489);
t555 = t406 * t601 * t610;
t607 = t376 * t588 + t604 * t555;
t396 = t420 + t405;
t598 = t396 * t405;
t381 = ((pkin(2) * t540 * t611 + t396 * t534 + t536 * t420) * t535 * t420 + (pkin(2) + t540) * t598) * t589;
t384 = ((-pkin(2) * t396 - t420 * t540) * t420 - pkin(2) * t598) * t589;
t480 = -rSges(2,1) * t514 + rSges(2,2) * t517;
t483 = rSges(2,1) * t517 + rSges(2,2) * t514;
t543 = (-t480 * t488 + t483 * t491) * pkin(1);
t438 = (t487 + 0.2e1 * t543) * m(2) + t553;
t441 = Icges(2,3) + (t503 + t543) * m(2);
t375 = -t381 * t441 - t384 * t438;
t602 = m(2) * (-t480 * t491 - t483 * t488);
t556 = t405 * t602 * t611;
t606 = t375 * t589 + t605 * t556;
t507 = t520 ^ 2;
t599 = m(3) * t507;
t595 = t459 * t605;
t594 = t460 * t605;
t593 = t461 * t604;
t592 = t462 * t604;
t591 = t463 * t603;
t590 = t464 * t603;
t586 = t605 * t535;
t585 = t604 * t535;
t584 = t603 * t535;
t486 = t503 * m(2) + Icges(2,3);
t583 = t486 * t535;
t581 = t537 * t507;
t580 = t420 ^ 2 * t602;
t579 = t421 ^ 2 * t601;
t578 = t422 ^ 2 * t600;
t577 = t444 * t586;
t576 = t445 * t585;
t575 = t446 * t584;
t574 = t447 * t586;
t573 = t448 * t585;
t572 = t449 * t584;
t571 = t605 * t583;
t569 = t604 * t583;
t567 = t603 * t583;
t552 = t580 * t586;
t551 = t579 * t585;
t550 = t578 * t584;
t525 = rSges(3,1);
t524 = rSges(3,2);
t425 = (t463 * t479 - t464 * t476) * t587;
t424 = (t461 * t478 - t462 * t475) * t588;
t423 = (t459 * t477 - t460 * t474) * t589;
t416 = (t446 * t479 - t449 * t476) * t566;
t415 = (t445 * t478 - t448 * t475) * t568;
t414 = (t444 * t477 - t447 * t474) * t570;
t413 = (t443 * t590 - t449 * t567) * t537;
t412 = (t443 * t591 - t446 * t567) * t537;
t411 = (t442 * t592 - t448 * t569) * t537;
t410 = (t442 * t593 - t445 * t569) * t537;
t409 = (t441 * t594 - t447 * t571) * t537;
t408 = (t441 * t595 - t444 * t571) * t537;
t404 = (t440 * t590 - t443 * t572) * t537;
t403 = (t440 * t591 - t443 * t575) * t537;
t402 = (t439 * t592 - t442 * t573) * t537;
t401 = (t439 * t593 - t442 * t576) * t537;
t400 = (t438 * t594 - t441 * t574) * t537;
t399 = (t438 * t595 - t441 * t577) * t537;
t392 = -t416 * t486 + t425 * t443;
t391 = -t415 * t486 + t424 * t442;
t390 = -t414 * t486 + t423 * t441;
t389 = -t416 * t443 + t425 * t440;
t388 = -t415 * t442 + t424 * t439;
t387 = -t414 * t441 + t423 * t438;
t380 = -t486 * t383 - t443 * t386;
t379 = -t486 * t382 - t442 * t385;
t378 = -t486 * t381 - t441 * t384;
t1 = [-t380 * t557 + t449 * t550 - t379 * t558 + t448 * t551 - t378 * t559 + t447 * t552 - (-t505 * t524 + t506 * t525) * t599 + t608 * t464 + t607 * t462 + t606 * t460 + (-(t404 * t590 - t413 * t572) * t479 - (t404 * t591 - t413 * t575) * t476 - (t402 * t592 - t411 * t573) * t478 - (t402 * t593 - t411 * t576) * t475 - (t400 * t594 - t409 * t574) * t477 - (t400 * t595 - t409 * t577) * t474) * t581; -t380 * t560 + t446 * t550 - t379 * t561 + t445 * t551 - t378 * t562 + t444 * t552 - (t505 * t525 + t506 * t524) * t599 + t608 * t463 + t607 * t461 + t606 * t459 + (-(t403 * t590 - t412 * t572) * t479 - (t403 * t591 - t412 * t575) * t476 - (t401 * t592 - t410 * t573) * t478 - (t401 * t593 - t410 * t576) * t475 - (t399 * t594 - t408 * t574) * t477 - (t399 * t595 - t408 * t577) * t474) * t581; t423 * t375 + t424 * t376 + t425 * t377 - t414 * t378 - t415 * t379 - t416 * t380 + (-(t389 * t590 - t392 * t572) * t479 - (t389 * t591 - t392 * t575) * t476 - (t388 * t592 - t391 * t573) * t478 - (t388 * t593 - t391 * t576) * t475 - (t387 * t594 - t390 * t574) * t477 - (t387 * t595 - t390 * t577) * t474) * t581 + (t414 * t580 + t415 * t579 + t416 * t578 + t423 * t556 + t424 * t555 + t425 * t554) * pkin(1);];
taucX  = t1;
