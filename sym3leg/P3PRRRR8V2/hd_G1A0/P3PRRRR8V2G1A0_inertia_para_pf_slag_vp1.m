% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G1A0
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:15
% EndTime: 2020-08-06 17:35:16
% DurationCPUTime: 1.25s
% Computational Cost: add. (3717->214), mult. (7107->408), div. (240->7), fcn. (6828->22), ass. (0->194)
t524 = cos(qJ(2,1));
t518 = sin(qJ(2,1));
t527 = pkin(7) + pkin(6);
t556 = t518 * t527;
t480 = pkin(2) * t524 + t556;
t506 = sin(pkin(8));
t508 = cos(pkin(8));
t493 = t527 * t524;
t477 = pkin(2) * t518 - t493;
t509 = cos(pkin(4));
t507 = sin(pkin(4));
t517 = sin(qJ(3,1));
t569 = t507 * t517;
t532 = pkin(3) * t569 - t477 * t509;
t607 = t480 * t508 + t506 * t532;
t522 = cos(qJ(2,2));
t516 = sin(qJ(2,2));
t558 = t516 * t527;
t479 = pkin(2) * t522 + t558;
t492 = t527 * t522;
t476 = pkin(2) * t516 - t492;
t515 = sin(qJ(3,2));
t571 = t507 * t515;
t533 = pkin(3) * t571 - t476 * t509;
t606 = t479 * t508 + t506 * t533;
t520 = cos(qJ(2,3));
t514 = sin(qJ(2,3));
t560 = t514 * t527;
t478 = pkin(2) * t520 + t560;
t491 = t527 * t520;
t475 = pkin(2) * t514 - t491;
t513 = sin(qJ(3,3));
t573 = t507 * t513;
t534 = pkin(3) * t573 - t475 * t509;
t605 = t478 * t508 + t506 * t534;
t603 = (m(3) * rSges(3,2));
t604 = -2 * rSges(3,1) * t603 + 2 * Icges(3,4);
t528 = pkin(2) * m(3);
t519 = cos(qJ(3,3));
t537 = rSges(3,1) * t519 - rSges(3,2) * t513;
t572 = t507 * t514;
t602 = m(3) * (t537 * t509 - (rSges(3,1) * t513 + rSges(3,2) * t519) * t572);
t521 = cos(qJ(3,2));
t536 = rSges(3,1) * t521 - rSges(3,2) * t515;
t570 = t507 * t516;
t601 = m(3) * (t536 * t509 - (rSges(3,1) * t515 + rSges(3,2) * t521) * t570);
t523 = cos(qJ(3,1));
t535 = rSges(3,1) * t523 - rSges(3,2) * t517;
t568 = t507 * t518;
t600 = m(3) * (t535 * t509 - (rSges(3,1) * t517 + rSges(3,2) * t523) * t568);
t503 = t519 ^ 2;
t599 = pkin(3) * t503;
t504 = t521 ^ 2;
t598 = pkin(3) * t504;
t505 = t523 ^ 2;
t597 = pkin(3) * t505;
t525 = pkin(6) + rSges(3,3);
t596 = t525 * m(3);
t510 = legFrame(3,3);
t495 = sin(t510);
t498 = cos(t510);
t457 = -t495 * t506 + t498 * t508;
t460 = t495 * t508 + t498 * t506;
t487 = pkin(3) * t519 + pkin(2);
t472 = t487 * t514 - t491;
t580 = (t487 * t520 + t560) * t509;
t430 = -t457 * t472 - t460 * t580;
t561 = t513 * t509;
t567 = t507 * t519;
t445 = 0.1e1 / (t472 * t567 + t487 * t561);
t595 = t430 * t445;
t511 = legFrame(2,3);
t496 = sin(t511);
t499 = cos(t511);
t458 = -t496 * t506 + t499 * t508;
t461 = t496 * t508 + t499 * t506;
t488 = pkin(3) * t521 + pkin(2);
t473 = t488 * t516 - t492;
t579 = (t488 * t522 + t558) * t509;
t431 = -t458 * t473 - t461 * t579;
t559 = t515 * t509;
t566 = t507 * t521;
t446 = 0.1e1 / (t473 * t566 + t488 * t559);
t594 = t431 * t446;
t512 = legFrame(1,3);
t497 = sin(t512);
t500 = cos(t512);
t459 = -t497 * t506 + t500 * t508;
t462 = t497 * t508 + t500 * t506;
t489 = pkin(3) * t523 + pkin(2);
t474 = t489 * t518 - t493;
t578 = (t489 * t524 + t556) * t509;
t432 = -t459 * t474 - t462 * t578;
t557 = t517 * t509;
t565 = t507 * t523;
t447 = 0.1e1 / (t474 * t565 + t489 * t557);
t593 = t432 * t447;
t433 = -t457 * t580 + t460 * t472;
t592 = t433 * t445;
t434 = -t458 * t579 + t461 * t473;
t591 = t434 * t446;
t435 = -t459 * t578 + t462 * t474;
t590 = t435 * t447;
t439 = 0.1e1 / (t572 * t599 + (pkin(3) * t561 + t475 * t507) * t519 + pkin(2) * t561);
t502 = m(1) + m(2) + m(3);
t589 = t439 * t502;
t440 = 0.1e1 / (t570 * t598 + (pkin(3) * t559 + t476 * t507) * t521 + pkin(2) * t559);
t588 = t440 * t502;
t441 = 0.1e1 / (t568 * t597 + (pkin(3) * t557 + t477 * t507) * t523 + pkin(2) * t557);
t587 = t441 * t502;
t531 = 0.1e1 / pkin(3);
t586 = t445 * t531;
t585 = t446 * t531;
t584 = t447 * t531;
t484 = m(2) * rSges(2,2) - t596;
t555 = m(2) * rSges(2,1) + t528;
t583 = ((m(3) * t537 + t555) * t520 - t484 * t514) * t507;
t582 = ((m(3) * t536 + t555) * t522 - t484 * t516) * t507;
t581 = ((m(3) * t535 + t555) * t524 - t484 * t518) * t507;
t529 = rSges(3,2) ^ 2;
t530 = rSges(3,1) ^ 2;
t574 = ((t529 + t530) * m(3) + Icges(3,3)) * t531;
t564 = t509 * t514;
t563 = t509 * t516;
t562 = t509 * t518;
t554 = pkin(2) * t573;
t553 = pkin(2) * t571;
t552 = pkin(2) * t569;
t551 = -0.2e1 * pkin(2) * t603;
t550 = t439 * t583;
t549 = t440 * t582;
t548 = t441 * t581;
t485 = -rSges(3,2) * t596 + Icges(3,6);
t486 = rSges(3,1) * t596 - Icges(3,5);
t454 = t485 * t519 - t486 * t513;
t547 = t454 * t586;
t546 = t445 * t574;
t455 = t485 * t521 - t486 * t515;
t545 = t455 * t585;
t544 = t446 * t574;
t456 = t485 * t523 - t486 * t517;
t543 = t456 * t584;
t542 = t447 * t574;
t442 = t506 * t478 - t508 * t534;
t463 = t506 * t564 - t508 * t520;
t466 = t506 * t520 + t508 * t564;
t406 = -(t463 * t498 + t466 * t495) * t599 + (-t442 * t495 + t498 * t605) * t519 + t460 * t554;
t424 = -t457 * t567 - (t457 * t564 + t460 * t520) * t513;
t541 = t586 * t602;
t388 = t406 * t589 + t424 * t550 + t433 * t541;
t443 = t506 * t479 - t508 * t533;
t464 = t506 * t563 - t508 * t522;
t467 = t506 * t522 + t508 * t563;
t407 = -(t464 * t499 + t467 * t496) * t598 + (-t443 * t496 + t499 * t606) * t521 + t461 * t553;
t425 = -t458 * t566 - (t458 * t563 + t461 * t522) * t515;
t540 = t585 * t601;
t389 = t407 * t588 + t425 * t549 + t434 * t540;
t444 = t506 * t480 - t508 * t532;
t465 = t506 * t562 - t508 * t524;
t468 = t506 * t524 + t508 * t562;
t408 = -(t465 * t500 + t468 * t497) * t597 + (-t444 * t497 + t500 * t607) * t523 + t462 * t552;
t426 = -t459 * t565 - (t459 * t562 + t462 * t524) * t517;
t539 = t584 * t600;
t390 = t408 * t587 + t426 * t548 + t435 * t539;
t409 = (-t463 * t495 + t466 * t498) * t599 + (t442 * t498 + t495 * t605) * t519 - t457 * t554;
t427 = -t460 * t567 - (-t457 * t520 + t460 * t564) * t513;
t391 = t409 * t589 + t427 * t550 + t430 * t541;
t410 = (-t464 * t496 + t467 * t499) * t598 + (t443 * t499 + t496 * t606) * t521 - t458 * t553;
t428 = -t461 * t566 - (-t458 * t522 + t461 * t563) * t515;
t392 = t410 * t588 + t428 * t549 + t431 * t540;
t411 = (-t465 * t497 + t468 * t500) * t597 + (t444 * t500 + t497 * t607) * t523 - t459 * t552;
t429 = -t462 * t565 - (-t459 * t524 + t462 * t562) * t517;
t393 = t411 * t587 + t429 * t548 + t432 * t539;
t538 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t525 ^ 2 + t529) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t501 = 0.2e1 * rSges(3,1) * t528;
t481 = (-t529 + t530) * m(3) + Icges(3,2) - Icges(3,1);
t438 = t481 * t505 + (t517 * t604 + t501) * t523 + t517 * t551 + t538;
t437 = t481 * t504 + (t515 * t604 + t501) * t521 + t515 * t551 + t538;
t436 = t481 * t503 + (t513 * t604 + t501) * t519 + t513 * t551 + t538;
t399 = t432 * t542 + (t411 * t600 + t429 * t456) * t441;
t398 = t431 * t544 + (t410 * t601 + t428 * t455) * t440;
t397 = t430 * t546 + (t409 * t602 + t427 * t454) * t439;
t396 = t435 * t542 + (t408 * t600 + t426 * t456) * t441;
t395 = t434 * t544 + (t407 * t601 + t425 * t455) * t440;
t394 = t433 * t546 + (t406 * t602 + t424 * t454) * t439;
t387 = t432 * t543 + (t411 * t581 + t429 * t438) * t441;
t386 = t431 * t545 + (t410 * t582 + t428 * t437) * t440;
t385 = t430 * t547 + (t409 * t583 + t427 * t436) * t439;
t384 = t435 * t543 + (t408 * t581 + t426 * t438) * t441;
t383 = t434 * t545 + (t407 * t582 + t425 * t437) * t440;
t382 = t433 * t547 + (t406 * t583 + t424 * t436) * t439;
t381 = t393 + t392 + t391;
t380 = t390 + t389 + t388;
t1 = [m(4) + (t384 * t426 + t390 * t408) * t441 + (t383 * t425 + t389 * t407) * t440 + (t382 * t424 + t388 * t406) * t439 + (t394 * t592 + t395 * t591 + t396 * t590) * t531, (t384 * t429 + t390 * t411) * t441 + (t383 * t428 + t389 * t410) * t440 + (t382 * t427 + t388 * t409) * t439 + (t394 * t595 + t395 * t594 + t396 * t593) * t531, t380; (t387 * t426 + t393 * t408) * t441 + (t386 * t425 + t392 * t407) * t440 + (t385 * t424 + t391 * t406) * t439 + (t397 * t592 + t398 * t591 + t399 * t590) * t531, m(4) + (t387 * t429 + t393 * t411) * t441 + (t386 * t428 + t392 * t410) * t440 + (t385 * t427 + t391 * t409) * t439 + (t397 * t595 + t398 * t594 + t399 * t593) * t531, t381; t380, t381, 0.3e1 * m(1) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
