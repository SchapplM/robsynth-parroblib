% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:08
% EndTime: 2020-08-06 17:49:09
% DurationCPUTime: 1.40s
% Computational Cost: add. (5502->250), mult. (11430->499), div. (432->4), fcn. (10116->22), ass. (0->200)
t612 = (m(3) * rSges(3,2));
t613 = -2 * rSges(3,1) * t612 + 2 * Icges(3,4);
t529 = pkin(2) * m(3);
t510 = cos(pkin(4));
t514 = sin(qJ(3,3));
t520 = cos(qJ(3,3));
t556 = rSges(3,1) * t520 - rSges(3,2) * t514;
t508 = sin(pkin(4));
t515 = sin(qJ(2,3));
t580 = t508 * t515;
t467 = t556 * t510 - (t514 * rSges(3,1) + t520 * rSges(3,2)) * t580;
t611 = m(3) * t467;
t516 = sin(qJ(3,2));
t522 = cos(qJ(3,2));
t555 = rSges(3,1) * t522 - rSges(3,2) * t516;
t517 = sin(qJ(2,2));
t578 = t508 * t517;
t468 = t555 * t510 - (t516 * rSges(3,1) + t522 * rSges(3,2)) * t578;
t610 = m(3) * t468;
t518 = sin(qJ(3,1));
t524 = cos(qJ(3,1));
t554 = rSges(3,1) * t524 - rSges(3,2) * t518;
t519 = sin(qJ(2,1));
t576 = t508 * t519;
t469 = t554 * t510 - (t518 * rSges(3,1) + t524 * rSges(3,2)) * t576;
t609 = m(3) * t469;
t532 = 0.1e1 / pkin(3);
t608 = m(3) * t532;
t607 = pkin(2) * t514;
t606 = pkin(2) * t516;
t605 = pkin(2) * t518;
t504 = t520 ^ 2;
t604 = pkin(3) * t504;
t505 = t522 ^ 2;
t603 = pkin(3) * t505;
t506 = t524 ^ 2;
t602 = pkin(3) * t506;
t601 = pkin(3) * t520;
t600 = pkin(3) * t522;
t599 = pkin(3) * t524;
t526 = pkin(6) + rSges(3,3);
t598 = t526 * m(3);
t521 = cos(qJ(2,3));
t528 = pkin(7) + pkin(6);
t482 = pkin(2) * t515 - t528 * t521;
t485 = pkin(2) * t521 + t515 * t528;
t507 = sin(pkin(8));
t509 = cos(pkin(8));
t568 = t510 * t521;
t575 = t509 * t510;
t446 = (t507 * t515 - t509 * t568) * t601 - t485 * t575 + t507 * t482;
t597 = t446 * t532;
t523 = cos(qJ(2,2));
t483 = pkin(2) * t517 - t528 * t523;
t486 = pkin(2) * t523 + t517 * t528;
t567 = t510 * t523;
t447 = (t507 * t517 - t509 * t567) * t600 - t486 * t575 + t507 * t483;
t596 = t447 * t532;
t525 = cos(qJ(2,1));
t484 = pkin(2) * t519 - t528 * t525;
t487 = pkin(2) * t525 + t519 * t528;
t566 = t510 * t525;
t448 = (t507 * t519 - t509 * t566) * t599 - t487 * t575 + t507 * t484;
t595 = t448 * t532;
t583 = t507 * t510;
t449 = (t507 * t568 + t509 * t515) * t601 + t485 * t583 + t482 * t509;
t594 = t449 * t532;
t450 = (t507 * t567 + t509 * t517) * t600 + t486 * t583 + t483 * t509;
t593 = t450 * t532;
t451 = (t507 * t566 + t509 * t519) * t599 + t487 * t583 + t484 * t509;
t592 = t451 * t532;
t491 = m(2) * rSges(2,2) - t598;
t562 = m(2) * rSges(2,1) + t529;
t591 = ((t556 * m(3) + t562) * t521 - t515 * t491) * t508;
t590 = ((t555 * m(3) + t562) * t523 - t517 * t491) * t508;
t589 = ((t554 * m(3) + t562) * t525 - t519 * t491) * t508;
t492 = -rSges(3,2) * t598 + Icges(3,6);
t493 = rSges(3,1) * t598 - Icges(3,5);
t473 = t492 * t520 - t493 * t514;
t588 = t473 * t532;
t474 = t492 * t522 - t493 * t516;
t587 = t474 * t532;
t475 = t492 * t524 - t493 * t518;
t586 = t475 * t532;
t530 = rSges(3,2) ^ 2;
t531 = rSges(3,1) ^ 2;
t585 = ((t530 + t531) * m(3) + Icges(3,3)) * t532;
t584 = t507 * t508;
t582 = t508 * t509;
t581 = t508 * t514;
t579 = t508 * t516;
t577 = t508 * t518;
t574 = t509 * t520;
t573 = t509 * t522;
t572 = t509 * t524;
t571 = t510 * t515;
t570 = t510 * t517;
t569 = t510 * t519;
t565 = t514 * t510;
t564 = t516 * t510;
t563 = t518 * t510;
t561 = t467 * t608;
t560 = t468 * t608;
t559 = t469 * t608;
t558 = -0.2e1 * pkin(2) * t612;
t557 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t526 ^ 2 + t530) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t553 = pkin(3) * t581 - t482 * t510;
t552 = pkin(3) * t579 - t483 * t510;
t551 = pkin(3) * t577 - t484 * t510;
t455 = -t507 * t485 + t553 * t509;
t470 = pkin(3) * t565 + t482 * t508;
t479 = t507 * t521 + t509 * t571;
t511 = legFrame(3,2);
t496 = sin(t511);
t499 = cos(t511);
t434 = -(t479 * t496 - t499 * t580) * t604 + (t455 * t496 + t499 * t470) * t520 + (t496 * t582 + t499 * t510) * t607;
t452 = 0.1e1 / (pkin(2) * t565 + t470 * t520 + t580 * t604);
t488 = (-t530 + t531) * m(3) + Icges(3,2) - Icges(3,1);
t502 = 0.2e1 * rSges(3,1) * t529;
t443 = t488 * t504 + (t514 * t613 + t502) * t520 + t514 * t558 + t557;
t476 = t507 * t571 - t509 * t521;
t461 = t514 * t476 + t520 * t584;
t541 = t443 * t461 + t449 * t588;
t407 = (t434 * t591 + t541 * t496) * t452;
t538 = t449 * t585 + t461 * t473;
t419 = (t434 * t611 + t538 * t496) * t452;
t550 = t407 * t461 + t419 * t594;
t435 = (t479 * t499 + t496 * t580) * t604 + (-t455 * t499 + t470 * t496) * t520 + (t510 * t496 - t499 * t582) * t607;
t408 = (t435 * t591 - t541 * t499) * t452;
t420 = (t435 * t611 - t538 * t499) * t452;
t549 = t408 * t461 + t420 * t594;
t456 = -t507 * t486 + t552 * t509;
t471 = pkin(3) * t564 + t483 * t508;
t480 = t507 * t523 + t509 * t570;
t512 = legFrame(2,2);
t497 = sin(t512);
t500 = cos(t512);
t436 = -(t480 * t497 - t500 * t578) * t603 + (t456 * t497 + t500 * t471) * t522 + (t497 * t582 + t500 * t510) * t606;
t453 = 0.1e1 / (pkin(2) * t564 + t471 * t522 + t578 * t603);
t444 = t488 * t505 + (t516 * t613 + t502) * t522 + t516 * t558 + t557;
t477 = t507 * t570 - t509 * t523;
t462 = t516 * t477 + t522 * t584;
t540 = t444 * t462 + t450 * t587;
t409 = (t436 * t590 + t540 * t497) * t453;
t537 = t450 * t585 + t462 * t474;
t421 = (t436 * t610 + t537 * t497) * t453;
t548 = t409 * t462 + t421 * t593;
t437 = (t480 * t500 + t497 * t578) * t603 + (-t456 * t500 + t471 * t497) * t522 + (t510 * t497 - t500 * t582) * t606;
t410 = (t437 * t590 - t540 * t500) * t453;
t422 = (t437 * t610 - t537 * t500) * t453;
t547 = t410 * t462 + t422 * t593;
t457 = -t507 * t487 + t551 * t509;
t472 = pkin(3) * t563 + t484 * t508;
t481 = t507 * t525 + t509 * t569;
t513 = legFrame(1,2);
t498 = sin(t513);
t501 = cos(t513);
t438 = -(t481 * t498 - t501 * t576) * t602 + (t457 * t498 + t501 * t472) * t524 + (t498 * t582 + t501 * t510) * t605;
t454 = 0.1e1 / (pkin(2) * t563 + t472 * t524 + t576 * t602);
t445 = t488 * t506 + (t518 * t613 + t502) * t524 + t518 * t558 + t557;
t478 = t507 * t569 - t509 * t525;
t463 = t518 * t478 + t524 * t584;
t539 = t445 * t463 + t451 * t586;
t411 = (t438 * t589 + t539 * t498) * t454;
t536 = t451 * t585 + t463 * t475;
t423 = (t438 * t609 + t536 * t498) * t454;
t546 = t411 * t463 + t423 * t592;
t439 = (t481 * t501 + t498 * t576) * t602 + (-t457 * t501 + t472 * t498) * t524 + (t510 * t498 - t501 * t582) * t605;
t412 = (t439 * t589 - t539 * t501) * t454;
t424 = (t439 * t609 - t536 * t501) * t454;
t545 = t412 * t463 + t424 * t592;
t440 = -t476 * t604 + t485 * t574 + (pkin(2) * t581 + t553 * t520) * t507;
t464 = -t514 * t479 - t508 * t574;
t425 = (t440 * t591 + t443 * t464 + t446 * t588) * t452;
t431 = (t440 * t611 + t446 * t585 + t464 * t473) * t452;
t544 = t425 * t461 + t431 * t594;
t441 = -t477 * t603 + t486 * t573 + (pkin(2) * t579 + t552 * t522) * t507;
t465 = -t516 * t480 - t508 * t573;
t426 = (t441 * t590 + t444 * t465 + t447 * t587) * t453;
t432 = (t441 * t610 + t447 * t585 + t465 * t474) * t453;
t543 = t426 * t462 + t432 * t593;
t442 = -t478 * t602 + t487 * t572 + (pkin(2) * t577 + t551 * t524) * t507;
t466 = -t518 * t481 - t508 * t572;
t427 = (t442 * t589 + t445 * t466 + t448 * t586) * t454;
t433 = (t442 * t609 + t448 * t585 + t466 * t475) * t454;
t542 = t427 * t463 + t433 * t592;
t535 = t449 * t561 + t461 * t591;
t534 = t450 * t560 + t462 * t590;
t533 = t451 * t559 + t463 * t589;
t503 = m(1) + m(2) + m(3);
t430 = (t442 * t503 + t448 * t559 + t466 * t589) * t454;
t429 = (t441 * t503 + t447 * t560 + t465 * t590) * t453;
t428 = (t440 * t503 + t446 * t561 + t464 * t591) * t452;
t418 = (t439 * t503 - t533 * t501) * t454;
t417 = (t438 * t503 + t533 * t498) * t454;
t416 = (t437 * t503 - t534 * t500) * t453;
t415 = (t436 * t503 + t534 * t497) * t453;
t414 = (t435 * t503 - t535 * t499) * t452;
t413 = (t434 * t503 + t535 * t496) * t452;
t1 = [m(4) + (t418 * t439 - t545 * t501) * t454 + (t416 * t437 - t547 * t500) * t453 + (t414 * t435 - t549 * t499) * t452, (t418 * t438 + t545 * t498) * t454 + (t416 * t436 + t547 * t497) * t453 + (t414 * t434 + t549 * t496) * t452, (t412 * t466 + t418 * t442 + t424 * t595) * t454 + (t410 * t465 + t416 * t441 + t422 * t596) * t453 + (t408 * t464 + t414 * t440 + t420 * t597) * t452; (t417 * t439 - t546 * t501) * t454 + (t415 * t437 - t548 * t500) * t453 + (t413 * t435 - t550 * t499) * t452, m(4) + (t417 * t438 + t546 * t498) * t454 + (t415 * t436 + t548 * t497) * t453 + (t413 * t434 + t550 * t496) * t452, (t411 * t466 + t417 * t442 + t423 * t595) * t454 + (t409 * t465 + t415 * t441 + t421 * t596) * t453 + (t407 * t464 + t413 * t440 + t419 * t597) * t452; (t430 * t439 - t542 * t501) * t454 + (t429 * t437 - t543 * t500) * t453 + (t428 * t435 - t544 * t499) * t452, (t430 * t438 + t542 * t498) * t454 + (t429 * t436 + t543 * t497) * t453 + (t428 * t434 + t544 * t496) * t452, m(4) + (t427 * t466 + t430 * t442 + t433 * t595) * t454 + (t426 * t465 + t429 * t441 + t432 * t596) * t453 + (t425 * t464 + t428 * t440 + t431 * t597) * t452;];
MX  = t1;
