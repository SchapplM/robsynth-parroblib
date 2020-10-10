% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G3A0
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:13
% EndTime: 2020-08-06 18:04:15
% DurationCPUTime: 1.39s
% Computational Cost: add. (5502->250), mult. (11430->501), div. (432->4), fcn. (10116->22), ass. (0->199)
t604 = (m(3) * rSges(3,2));
t605 = -2 * rSges(3,1) * t604 + 2 * Icges(3,4);
t522 = pkin(2) * m(3);
t503 = cos(pkin(4));
t507 = sin(qJ(3,3));
t513 = cos(qJ(3,3));
t549 = rSges(3,1) * t513 - rSges(3,2) * t507;
t501 = sin(pkin(4));
t508 = sin(qJ(2,3));
t573 = t501 * t508;
t460 = t549 * t503 - (t507 * rSges(3,1) + t513 * rSges(3,2)) * t573;
t603 = m(3) * t460;
t509 = sin(qJ(3,2));
t515 = cos(qJ(3,2));
t548 = rSges(3,1) * t515 - rSges(3,2) * t509;
t510 = sin(qJ(2,2));
t571 = t501 * t510;
t461 = t548 * t503 - (t509 * rSges(3,1) + t515 * rSges(3,2)) * t571;
t602 = m(3) * t461;
t511 = sin(qJ(3,1));
t517 = cos(qJ(3,1));
t547 = rSges(3,1) * t517 - rSges(3,2) * t511;
t512 = sin(qJ(2,1));
t569 = t501 * t512;
t462 = t547 * t503 - (t511 * rSges(3,1) + t517 * rSges(3,2)) * t569;
t601 = m(3) * t462;
t525 = 0.1e1 / pkin(3);
t600 = m(3) * t525;
t599 = pkin(2) * t507;
t598 = pkin(2) * t509;
t597 = pkin(2) * t511;
t497 = t513 ^ 2;
t596 = pkin(3) * t497;
t498 = t515 ^ 2;
t595 = pkin(3) * t498;
t499 = t517 ^ 2;
t594 = pkin(3) * t499;
t593 = pkin(3) * t513;
t592 = pkin(3) * t515;
t591 = pkin(3) * t517;
t519 = pkin(6) + rSges(3,3);
t590 = t519 * m(3);
t514 = cos(qJ(2,3));
t521 = pkin(7) + pkin(6);
t475 = pkin(2) * t508 - t521 * t514;
t478 = pkin(2) * t514 + t508 * t521;
t500 = sin(pkin(8));
t502 = cos(pkin(8));
t561 = t503 * t514;
t565 = t502 * t503;
t439 = (t500 * t508 - t502 * t561) * t593 - t478 * t565 + t475 * t500;
t589 = t439 * t525;
t516 = cos(qJ(2,2));
t476 = pkin(2) * t510 - t521 * t516;
t479 = pkin(2) * t516 + t510 * t521;
t560 = t503 * t516;
t440 = (t500 * t510 - t502 * t560) * t592 - t479 * t565 + t476 * t500;
t588 = t440 * t525;
t518 = cos(qJ(2,1));
t477 = pkin(2) * t512 - t521 * t518;
t480 = pkin(2) * t518 + t512 * t521;
t559 = t503 * t518;
t441 = (t500 * t512 - t502 * t559) * t591 - t480 * t565 + t477 * t500;
t587 = t441 * t525;
t575 = t500 * t503;
t442 = (t500 * t561 + t502 * t508) * t593 + t478 * t575 + t475 * t502;
t586 = t442 * t525;
t443 = (t500 * t560 + t502 * t510) * t592 + t479 * t575 + t476 * t502;
t585 = t443 * t525;
t444 = (t500 * t559 + t502 * t512) * t591 + t480 * t575 + t477 * t502;
t584 = t444 * t525;
t484 = m(2) * rSges(2,2) - t590;
t555 = m(2) * rSges(2,1) + t522;
t583 = ((t549 * m(3) + t555) * t514 - t508 * t484) * t501;
t582 = ((t548 * m(3) + t555) * t516 - t510 * t484) * t501;
t581 = ((t547 * m(3) + t555) * t518 - t512 * t484) * t501;
t485 = -rSges(3,2) * t590 + Icges(3,6);
t486 = rSges(3,1) * t590 - Icges(3,5);
t466 = t485 * t513 - t486 * t507;
t580 = t466 * t525;
t467 = t485 * t515 - t486 * t509;
t579 = t467 * t525;
t468 = t485 * t517 - t486 * t511;
t578 = t468 * t525;
t523 = rSges(3,2) ^ 2;
t524 = rSges(3,1) ^ 2;
t577 = ((t523 + t524) * m(3) + Icges(3,3)) * t525;
t576 = t500 * t501;
t574 = t501 * t507;
t572 = t501 * t509;
t570 = t501 * t511;
t568 = t501 * t513;
t567 = t501 * t515;
t566 = t501 * t517;
t564 = t503 * t508;
t563 = t503 * t510;
t562 = t503 * t512;
t558 = t507 * t503;
t557 = t509 * t503;
t556 = t511 * t503;
t554 = t460 * t600;
t553 = t461 * t600;
t552 = t462 * t600;
t551 = -0.2e1 * pkin(2) * t604;
t550 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t519 ^ 2 + t523) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t546 = pkin(3) * t574 - t475 * t503;
t545 = pkin(3) * t572 - t476 * t503;
t544 = pkin(3) * t570 - t477 * t503;
t448 = t502 * t478 + t546 * t500;
t463 = pkin(3) * t558 + t501 * t475;
t469 = t500 * t564 - t502 * t514;
t504 = legFrame(3,2);
t489 = sin(t504);
t492 = cos(t504);
t427 = -(t469 * t492 - t489 * t573) * t596 + (t448 * t492 + t489 * t463) * t513 + (t503 * t489 + t492 * t576) * t599;
t445 = 0.1e1 / (pkin(2) * t558 + t463 * t513 + t573 * t596);
t481 = (-t523 + t524) * m(3) + Icges(3,2) - Icges(3,1);
t495 = 0.2e1 * rSges(3,1) * t522;
t436 = t481 * t497 + (t507 * t605 + t495) * t513 + t507 * t551 + t550;
t472 = t500 * t514 + t502 * t564;
t457 = t507 * t472 + t502 * t568;
t534 = -t436 * t457 + t439 * t580;
t400 = (t427 * t583 + t534 * t492) * t445;
t531 = t439 * t577 - t457 * t466;
t412 = (t427 * t603 + t531 * t492) * t445;
t543 = -t400 * t457 + t412 * t589;
t449 = t502 * t479 + t545 * t500;
t464 = pkin(3) * t557 + t501 * t476;
t470 = t500 * t563 - t502 * t516;
t505 = legFrame(2,2);
t490 = sin(t505);
t493 = cos(t505);
t428 = -(t470 * t493 - t490 * t571) * t595 + (t449 * t493 + t490 * t464) * t515 + (t503 * t490 + t493 * t576) * t598;
t446 = 0.1e1 / (pkin(2) * t557 + t464 * t515 + t571 * t595);
t437 = t481 * t498 + (t509 * t605 + t495) * t515 + t509 * t551 + t550;
t473 = t500 * t516 + t502 * t563;
t458 = t509 * t473 + t502 * t567;
t533 = -t437 * t458 + t440 * t579;
t401 = (t428 * t582 + t533 * t493) * t446;
t530 = t440 * t577 - t458 * t467;
t413 = (t428 * t602 + t530 * t493) * t446;
t542 = -t401 * t458 + t413 * t588;
t450 = t502 * t480 + t544 * t500;
t465 = pkin(3) * t556 + t501 * t477;
t471 = t500 * t562 - t502 * t518;
t506 = legFrame(1,2);
t491 = sin(t506);
t494 = cos(t506);
t429 = -(t471 * t494 - t491 * t569) * t594 + (t450 * t494 + t491 * t465) * t517 + (t503 * t491 + t494 * t576) * t597;
t447 = 0.1e1 / (pkin(2) * t556 + t465 * t517 + t569 * t594);
t438 = t481 * t499 + (t511 * t605 + t495) * t517 + t511 * t551 + t550;
t474 = t500 * t518 + t502 * t562;
t459 = t511 * t474 + t502 * t566;
t532 = -t438 * t459 + t441 * t578;
t402 = (t429 * t581 + t532 * t494) * t447;
t529 = t441 * t577 - t459 * t468;
t414 = (t429 * t601 + t529 * t494) * t447;
t541 = -t402 * t459 + t414 * t587;
t430 = (t469 * t489 + t492 * t573) * t596 + (-t448 * t489 + t492 * t463) * t513 + (-t489 * t576 + t492 * t503) * t599;
t403 = (t430 * t583 - t534 * t489) * t445;
t415 = (t430 * t603 - t531 * t489) * t445;
t540 = -t403 * t457 + t415 * t589;
t431 = (t470 * t490 + t493 * t571) * t595 + (-t449 * t490 + t493 * t464) * t515 + (-t490 * t576 + t493 * t503) * t598;
t404 = (t431 * t582 - t533 * t490) * t446;
t416 = (t431 * t602 - t530 * t490) * t446;
t539 = -t404 * t458 + t416 * t588;
t432 = (t471 * t491 + t494 * t569) * t594 + (-t450 * t491 + t494 * t465) * t517 + (-t491 * t576 + t494 * t503) * t597;
t405 = (t432 * t581 - t532 * t491) * t447;
t417 = (t432 * t601 - t529 * t491) * t447;
t538 = -t405 * t459 + t417 * t587;
t433 = -t472 * t596 - t478 * t500 * t513 + (pkin(2) * t574 + t546 * t513) * t502;
t454 = t507 * t469 + t500 * t568;
t418 = (t433 * t583 + t436 * t454 + t442 * t580) * t445;
t424 = (t433 * t603 + t442 * t577 + t454 * t466) * t445;
t537 = -t418 * t457 + t424 * t589;
t434 = -t473 * t595 - t479 * t500 * t515 + (pkin(2) * t572 + t545 * t515) * t502;
t455 = t509 * t470 + t500 * t567;
t419 = (t434 * t582 + t437 * t455 + t443 * t579) * t446;
t425 = (t434 * t602 + t443 * t577 + t455 * t467) * t446;
t536 = -t419 * t458 + t425 * t588;
t435 = -t474 * t594 - t480 * t500 * t517 + (pkin(2) * t570 + t544 * t517) * t502;
t456 = t511 * t471 + t500 * t566;
t420 = (t435 * t581 + t438 * t456 + t444 * t578) * t447;
t426 = (t435 * t601 + t444 * t577 + t456 * t468) * t447;
t535 = -t420 * t459 + t426 * t587;
t528 = t439 * t554 - t457 * t583;
t527 = t440 * t553 - t458 * t582;
t526 = t441 * t552 - t459 * t581;
t496 = m(1) + m(2) + m(3);
t423 = (t435 * t496 + t444 * t552 + t456 * t581) * t447;
t422 = (t434 * t496 + t443 * t553 + t455 * t582) * t446;
t421 = (t433 * t496 + t442 * t554 + t454 * t583) * t445;
t411 = (t432 * t496 - t526 * t491) * t447;
t410 = (t431 * t496 - t527 * t490) * t446;
t409 = (t430 * t496 - t528 * t489) * t445;
t408 = (t429 * t496 + t526 * t494) * t447;
t407 = (t428 * t496 + t527 * t493) * t446;
t406 = (t427 * t496 + t528 * t492) * t445;
t1 = [m(4) + (t408 * t429 + t541 * t494) * t447 + (t407 * t428 + t542 * t493) * t446 + (t406 * t427 + t543 * t492) * t445, (t408 * t432 - t541 * t491) * t447 + (t407 * t431 - t542 * t490) * t446 + (t406 * t430 - t543 * t489) * t445, (t402 * t456 + t408 * t435 + t414 * t584) * t447 + (t401 * t455 + t407 * t434 + t413 * t585) * t446 + (t400 * t454 + t406 * t433 + t412 * t586) * t445; (t411 * t429 + t538 * t494) * t447 + (t410 * t428 + t539 * t493) * t446 + (t409 * t427 + t540 * t492) * t445, m(4) + (t411 * t432 - t538 * t491) * t447 + (t410 * t431 - t539 * t490) * t446 + (t409 * t430 - t540 * t489) * t445, (t405 * t456 + t411 * t435 + t417 * t584) * t447 + (t404 * t455 + t410 * t434 + t416 * t585) * t446 + (t403 * t454 + t409 * t433 + t415 * t586) * t445; (t423 * t429 + t535 * t494) * t447 + (t422 * t428 + t536 * t493) * t446 + (t421 * t427 + t537 * t492) * t445, (t423 * t432 - t535 * t491) * t447 + (t422 * t431 - t536 * t490) * t446 + (t421 * t430 - t537 * t489) * t445, m(4) + (t420 * t456 + t423 * t435 + t426 * t584) * t447 + (t419 * t455 + t422 * t434 + t425 * t585) * t446 + (t418 * t454 + t421 * t433 + t424 * t586) * t445;];
MX  = t1;
