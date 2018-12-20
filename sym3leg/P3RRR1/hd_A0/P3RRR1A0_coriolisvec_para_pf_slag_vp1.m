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
%   mass of all robot links (including platform)
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
% Datum: 2018-12-20 18:13
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-12-20 18:13:07
% EndTime: 2018-12-20 18:13:09
% DurationCPUTime: 2.21s
% Computational Cost: add. (10090->202), mult. (11590->404), div. (2436->5), fcn. (11938->20), ass. (0->199)
t537 = 0.1e1 / pkin(1);
t509 = qJ(1,2) + qJ(2,2);
t489 = sin(t509);
t492 = cos(t509);
t515 = sin(qJ(1,2));
t518 = cos(qJ(1,2));
t607 = 0.1e1 / (t489 * t518 - t515 * t492);
t591 = t607 * t537;
t508 = qJ(1,3) + qJ(2,3);
t488 = sin(t508);
t491 = cos(t508);
t514 = sin(qJ(1,3));
t517 = cos(qJ(1,3));
t608 = 0.1e1 / (t488 * t517 - t514 * t491);
t592 = t608 * t537;
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
t585 = t535 * t537;
t567 = t608 * t585;
t559 = t447 * t567;
t459 = t497 * t488 + t494 * t491;
t444 = pkin(1) * (t494 * t517 + t497 * t514) + t459 * pkin(2);
t562 = t444 * t567;
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
t565 = t607 * t585;
t558 = t448 * t565;
t461 = t498 * t489 + t495 * t492;
t445 = pkin(1) * (t495 * t518 + t498 * t515) + t461 * pkin(2);
t561 = t445 * t565;
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
t606 = 0.1e1 / (t490 * t519 - t516 * t493);
t563 = t606 * t585;
t557 = t449 * t563;
t463 = t499 * t490 + t496 * t493;
t446 = pkin(1) * (t496 * t519 + t499 * t516) + t463 * pkin(2);
t560 = t446 * t563;
t407 = -t455 * t560 - t458 * t557;
t573 = t460 * t592;
t574 = t459 * t592;
t420 = t453 * t574 + t456 * t573;
t611 = 0.2e1 * t420 + t405;
t571 = t462 * t591;
t572 = t461 * t591;
t421 = t454 * t572 + t457 * t571;
t610 = 0.2e1 * t421 + t406;
t590 = t606 * t537;
t569 = t464 * t590;
t570 = t463 * t590;
t422 = t455 * t570 + t458 * t569;
t609 = 0.2e1 * t422 + t407;
t538 = (t490 * t516 + t493 * t519) * pkin(1);
t539 = (t489 * t515 + t492 * t518) * pkin(1);
t540 = (t488 * t514 + t491 * t517) * pkin(1);
t507 = t520 ^ 2;
t480 = -rSges(2,1) * t514 + rSges(2,2) * t517;
t483 = rSges(2,1) * t517 + rSges(2,2) * t514;
t605 = m(2) * (-t480 * t491 - t483 * t488);
t481 = -rSges(2,1) * t515 + rSges(2,2) * t518;
t484 = rSges(2,1) * t518 + rSges(2,2) * t515;
t604 = m(2) * (-t481 * t492 - t484 * t489);
t482 = -rSges(2,1) * t516 + rSges(2,2) * t519;
t485 = rSges(2,1) * t519 + rSges(2,2) * t516;
t603 = m(2) * (-t482 * t493 - t485 * t490);
t602 = m(3) * t507;
t396 = t420 + t405;
t601 = t396 * t405;
t397 = t421 + t406;
t600 = t397 * t406;
t398 = t422 + t407;
t599 = t398 * t407;
t598 = t459 * t608;
t597 = t460 * t608;
t596 = t461 * t607;
t595 = t462 * t607;
t594 = t463 * t606;
t593 = t464 * t606;
t589 = t608 * t535;
t588 = t607 * t535;
t587 = t606 * t535;
t503 = rSges(2,1) ^ 2 + rSges(2,2) ^ 2;
t486 = t503 * m(2) + Icges(2,3);
t586 = t486 * t535;
t584 = t537 * t507;
t583 = t420 ^ 2 * t605;
t582 = t421 ^ 2 * t604;
t581 = t422 ^ 2 * t603;
t580 = t444 * t589;
t579 = t445 * t588;
t578 = t446 * t587;
t577 = t447 * t589;
t576 = t448 * t588;
t575 = t449 * t587;
t568 = t608 * t586;
t566 = t607 * t586;
t564 = t606 * t586;
t556 = t405 * t605 * t611;
t555 = t406 * t604 * t610;
t554 = t407 * t603 * t609;
t553 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(1,3) + Icges(2,3);
t552 = t583 * t589;
t551 = t582 * t588;
t550 = t581 * t587;
t546 = t608 * t556;
t545 = t607 * t555;
t544 = t606 * t554;
t543 = (-t480 * t488 + t483 * t491) * pkin(1);
t542 = (-t481 * t489 + t484 * t492) * pkin(1);
t541 = (-t482 * t490 + t485 * t493) * pkin(1);
t536 = pkin(1) ^ 2;
t534 = pkin(2) ^ 2;
t525 = rSges(3,1);
t524 = rSges(3,2);
t487 = t536 + t503;
t443 = Icges(2,3) + (t503 + t541) * m(2);
t442 = Icges(2,3) + (t503 + t542) * m(2);
t441 = Icges(2,3) + (t503 + t543) * m(2);
t440 = (t487 + 0.2e1 * t541) * m(2) + t553;
t439 = (t487 + 0.2e1 * t542) * m(2) + t553;
t438 = (t487 + 0.2e1 * t543) * m(2) + t553;
t425 = (t463 * t479 - t464 * t476) * t590;
t424 = (t461 * t478 - t462 * t475) * t591;
t423 = (t459 * t477 - t460 * t474) * t592;
t416 = (t446 * t479 - t449 * t476) * t563;
t415 = (t445 * t478 - t448 * t475) * t565;
t414 = (t444 * t477 - t447 * t474) * t567;
t413 = (t443 * t593 - t449 * t564) * t537;
t412 = (t443 * t594 - t446 * t564) * t537;
t411 = (t442 * t595 - t448 * t566) * t537;
t410 = (t442 * t596 - t445 * t566) * t537;
t409 = (t441 * t597 - t447 * t568) * t537;
t408 = (t441 * t598 - t444 * t568) * t537;
t404 = (t440 * t593 - t443 * t575) * t537;
t403 = (t440 * t594 - t443 * t578) * t537;
t402 = (t439 * t595 - t442 * t576) * t537;
t401 = (t439 * t596 - t442 * t579) * t537;
t400 = (t438 * t597 - t441 * t577) * t537;
t399 = (t438 * t598 - t441 * t580) * t537;
t392 = -t416 * t486 + t425 * t443;
t391 = -t415 * t486 + t424 * t442;
t390 = -t414 * t486 + t423 * t441;
t389 = -t416 * t443 + t425 * t440;
t388 = -t415 * t442 + t424 * t439;
t387 = -t414 * t441 + t423 * t438;
t386 = (-pkin(2) * t599 + (-t398 * pkin(2) - t422 * t538) * t422) * t590;
t385 = ((-pkin(2) * t397 - t421 * t539) * t421 - pkin(2) * t600) * t591;
t384 = ((-pkin(2) * t396 - t420 * t540) * t420 - pkin(2) * t601) * t592;
t383 = ((pkin(2) * t538 * t609 + t398 * t534 + t536 * t422) * t535 * t422 + (pkin(2) + t538) * t599) * t590;
t382 = ((pkin(2) * t539 * t610 + t397 * t534 + t536 * t421) * t535 * t421 + (pkin(2) + t539) * t600) * t591;
t381 = ((pkin(2) * t540 * t611 + t396 * t534 + t536 * t420) * t535 * t420 + (pkin(2) + t540) * t601) * t592;
t380 = -t486 * t383 - t443 * t386;
t379 = -t486 * t382 - t442 * t385;
t378 = -t486 * t381 - t441 * t384;
t377 = -t443 * t383 - t440 * t386;
t376 = -t442 * t382 - t439 * t385;
t375 = -t441 * t381 - t438 * t384;
t1 = [t377 * t569 - t380 * t557 + t464 * t544 + t449 * t550 + t376 * t571 - t379 * t558 + t462 * t545 + t448 * t551 + t375 * t573 - t378 * t559 + t460 * t546 + t447 * t552 - (-t505 * t524 + t506 * t525) * t602 + (-(t404 * t593 - t413 * t575) * t479 - (t404 * t594 - t413 * t578) * t476 - (t402 * t595 - t411 * t576) * t478 - (t402 * t596 - t411 * t579) * t475 - (t400 * t597 - t409 * t577) * t477 - (t400 * t598 - t409 * t580) * t474) * t584; t377 * t570 - t380 * t560 + t463 * t544 + t446 * t550 + t376 * t572 - t379 * t561 + t461 * t545 + t445 * t551 + t375 * t574 - t378 * t562 + t459 * t546 + t444 * t552 - (t505 * t525 + t506 * t524) * t602 + (-(t403 * t593 - t412 * t575) * t479 - (t403 * t594 - t412 * t578) * t476 - (t401 * t595 - t410 * t576) * t478 - (t401 * t596 - t410 * t579) * t475 - (t399 * t597 - t408 * t577) * t477 - (t399 * t598 - t408 * t580) * t474) * t584; t423 * t375 + t424 * t376 + t425 * t377 - t414 * t378 - t415 * t379 - t416 * t380 + (-(t389 * t593 - t392 * t575) * t479 - (t389 * t594 - t392 * t578) * t476 - (t388 * t595 - t391 * t576) * t478 - (t388 * t596 - t391 * t579) * t475 - (t387 * t597 - t390 * t577) * t477 - (t387 * t598 - t390 * t580) * t474) * t584 + (t414 * t583 + t415 * t582 + t416 * t581 + t423 * t556 + t424 * t555 + t425 * t554) * pkin(1);];
taucX  = t1;
