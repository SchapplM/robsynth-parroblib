% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:14:40
% EndTime: 2020-08-06 19:14:42
% DurationCPUTime: 1.80s
% Computational Cost: add. (4293->277), mult. (5064->458), div. (396->6), fcn. (3423->18), ass. (0->195)
t641 = 2 * rSges(3,3);
t644 = (t641 + qJ(3,3)) * qJ(3,3);
t643 = (t641 + qJ(3,2)) * qJ(3,2);
t642 = (t641 + qJ(3,1)) * qJ(3,1);
t564 = pkin(2) + pkin(3);
t640 = m(2) * rSges(2,1);
t562 = m(2) * rSges(2,2);
t561 = rSges(2,3) + pkin(5);
t560 = pkin(2) + rSges(3,1);
t547 = sin(qJ(2,3));
t613 = t547 * qJ(3,3);
t504 = pkin(1) + t613;
t553 = cos(qJ(2,3));
t604 = t564 * t553;
t492 = 0.1e1 / (t504 + t604);
t639 = m(3) * t492;
t549 = sin(qJ(2,2));
t611 = t549 * qJ(3,2);
t506 = pkin(1) + t611;
t555 = cos(qJ(2,2));
t603 = t564 * t555;
t493 = 0.1e1 / (t506 + t603);
t638 = m(3) * t493;
t551 = sin(qJ(2,1));
t609 = t551 * qJ(3,1);
t508 = pkin(1) + t609;
t557 = cos(qJ(2,1));
t602 = t564 * t557;
t494 = 0.1e1 / (t508 + t602);
t637 = m(3) * t494;
t559 = pkin(5) + rSges(3,2);
t636 = m(3) * t559;
t635 = pkin(1) * t547;
t634 = pkin(1) * t549;
t633 = pkin(1) * t551;
t541 = -qJ(3,3) - rSges(3,3);
t632 = t541 * m(3);
t542 = -qJ(3,2) - rSges(3,3);
t631 = t542 * m(3);
t543 = -qJ(3,1) - rSges(3,3);
t630 = t543 * m(3);
t629 = t560 * m(3);
t544 = legFrame(3,3);
t521 = sin(t544);
t524 = cos(t544);
t548 = sin(qJ(1,3));
t554 = cos(qJ(1,3));
t480 = t521 * t548 - t524 * t554;
t563 = pkin(5) - pkin(6);
t518 = t563 * t554;
t487 = t548 * t504 - t518;
t512 = t548 * t563;
t589 = t504 * t554 + t512;
t458 = t480 * t604 + t487 * t521 - t589 * t524;
t628 = t458 * t553;
t481 = t521 * t554 + t524 * t548;
t459 = t481 * t604 + t487 * t524 + t521 * t589;
t627 = t459 * t553;
t545 = legFrame(2,3);
t522 = sin(t545);
t525 = cos(t545);
t550 = sin(qJ(1,2));
t556 = cos(qJ(1,2));
t482 = t522 * t550 - t525 * t556;
t519 = t563 * t556;
t489 = t550 * t506 - t519;
t513 = t550 * t563;
t587 = t506 * t556 + t513;
t460 = t482 * t603 + t489 * t522 - t587 * t525;
t626 = t460 * t555;
t483 = t522 * t556 + t525 * t550;
t461 = t483 * t603 + t489 * t525 + t522 * t587;
t625 = t461 * t555;
t546 = legFrame(1,3);
t523 = sin(t546);
t526 = cos(t546);
t552 = sin(qJ(1,1));
t558 = cos(qJ(1,1));
t484 = t523 * t552 - t526 * t558;
t520 = t563 * t558;
t491 = t552 * t508 - t520;
t514 = t552 * t563;
t585 = t508 * t558 + t514;
t462 = t484 * t602 + t491 * t523 - t585 * t526;
t624 = t462 * t557;
t485 = t523 * t558 + t526 * t552;
t463 = t485 * t602 + t491 * t526 + t523 * t585;
t623 = t463 * t557;
t571 = rSges(3,3) ^ 2;
t584 = pkin(2) ^ 2 + t571 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t572 = rSges(2,2) ^ 2;
t574 = rSges(2,1) ^ 2;
t592 = Icges(3,2) + Icges(2,3) + (t572 + t574) * m(2);
t470 = (t584 + t644) * m(3) + t592;
t495 = -t553 * qJ(3,3) + t564 * t547;
t566 = 1 / qJ(3,3);
t464 = (t470 * t547 - t495 * t629) * t566;
t622 = t464 * t553;
t471 = (t584 + t643) * m(3) + t592;
t496 = -t555 * qJ(3,2) + t564 * t549;
t568 = 1 / qJ(3,2);
t465 = (t471 * t549 - t496 * t629) * t568;
t621 = t465 * t555;
t472 = (t584 + t642) * m(3) + t592;
t497 = -t557 * qJ(3,1) + t564 * t551;
t570 = 1 / qJ(3,1);
t466 = (t472 * t551 - t497 * t629) * t570;
t620 = t466 * t557;
t479 = t559 * t629 + t561 * t640 - Icges(3,4) - Icges(2,5);
t577 = -t561 * t562 + Icges(2,6) - Icges(3,6);
t467 = (-t559 * t632 + t577) * t553 - t547 * t479;
t619 = t467 * t553;
t468 = (-t559 * t631 + t577) * t555 - t549 * t479;
t618 = t468 * t555;
t469 = (-t559 * t630 + t577) * t557 - t551 * t479;
t617 = t469 * t557;
t616 = t470 * t553;
t615 = t471 * t555;
t614 = t472 * t557;
t612 = t547 * t559;
t610 = t549 * t559;
t608 = t551 * t559;
t607 = t553 * t560;
t606 = t555 * t560;
t605 = t557 * t560;
t601 = -Icges(2,1) - Icges(3,1);
t600 = rSges(3,3) - t560;
t599 = rSges(3,3) + t560;
t598 = m(3) * t612;
t597 = m(3) * t610;
t596 = m(3) * t608;
t537 = t553 ^ 2;
t595 = (qJ(3,3) + t564) * (-qJ(3,3) + t564) * t537;
t538 = t555 ^ 2;
t594 = (qJ(3,2) + t564) * (-qJ(3,2) + t564) * t538;
t539 = t557 ^ 2;
t593 = (qJ(3,1) + t564) * (-qJ(3,1) + t564) * t539;
t576 = pkin(1) ^ 2;
t591 = t559 ^ 2 + t571 + t576;
t503 = pkin(1) + 0.2e1 * t613;
t590 = t503 * t554 + t512;
t505 = pkin(1) + 0.2e1 * t611;
t588 = t505 * t556 + t513;
t507 = pkin(1) + 0.2e1 * t609;
t586 = t507 * t558 + t514;
t583 = Icges(1,3) + (t561 ^ 2 + t572 + t576) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t601;
t582 = (-t572 + t574) * m(2) + Icges(2,2) + Icges(3,3) + t601;
t581 = -rSges(2,1) * t562 + Icges(2,4) - Icges(3,5);
t509 = qJ(3,3) + t635;
t580 = t509 * t554 + t547 * t512;
t510 = qJ(3,2) + t634;
t579 = t510 * t556 + t549 * t513;
t511 = qJ(3,1) + t633;
t578 = t511 * t558 + t551 * t514;
t498 = 0.2e1 * (t629 + t640) * pkin(1);
t490 = t552 * t507 - t520;
t488 = t550 * t505 - t519;
t486 = t548 * t503 - t518;
t478 = t552 * t511 - t551 * t520;
t477 = t550 * t510 - t549 * t519;
t476 = t548 * t509 - t547 * t518;
t475 = (-t551 * t560 + t497) * t570 * m(3);
t474 = (-t549 * t560 + t496) * t568 * m(3);
t473 = (-t547 * t560 + t495) * t566 * m(3);
t457 = (t497 * t636 + t469) * t570 * t551;
t456 = (t496 * t636 + t468) * t568 * t549;
t455 = (t495 * t636 + t467) * t566 * t547;
t454 = (-(qJ(3,1) + t599) * (qJ(3,1) + t600) * m(3) + t582) * t539 + (0.2e1 * (-t543 * t629 + t581) * t551 + t498) * t557 - 0.2e1 * (t562 + t630) * t633 + (t591 + t642) * m(3) + t583;
t453 = (-(qJ(3,2) + t599) * (qJ(3,2) + t600) * m(3) + t582) * t538 + (0.2e1 * (-t542 * t629 + t581) * t549 + t498) * t555 - 0.2e1 * (t562 + t631) * t634 + (t591 + t643) * m(3) + t583;
t452 = (-(qJ(3,3) + t599) * (qJ(3,3) + t600) * m(3) + t582) * t537 + (0.2e1 * (-t541 * t629 + t581) * t547 + t498) * t553 - 0.2e1 * (t562 + t632) * t635 + (t591 + t644) * m(3) + t583;
t451 = t485 * t593 + (t490 * t526 + t586 * t523) * t602 + (t478 * t526 + t578 * t523) * qJ(3,1);
t450 = -t484 * t593 - (t523 * t490 - t586 * t526) * t602 - (t523 * t478 - t578 * t526) * qJ(3,1);
t449 = t483 * t594 + (t488 * t525 + t588 * t522) * t603 + (t477 * t525 + t579 * t522) * qJ(3,2);
t448 = -t482 * t594 - (t522 * t488 - t588 * t525) * t603 - (t522 * t477 - t579 * t525) * qJ(3,2);
t447 = t481 * t595 + (t486 * t524 + t590 * t521) * t604 + (t476 * t524 + t580 * t521) * qJ(3,3);
t446 = -t480 * t595 - (t521 * t486 - t590 * t524) * t604 - (t521 * t476 - t580 * t524) * qJ(3,3);
t445 = (-t485 * t608 + (t462 * t605 + t450) * t570) * t637;
t444 = (-t484 * t608 + (-t463 * t605 + t451) * t570) * t637;
t443 = (-t483 * t610 + (t460 * t606 + t448) * t568) * t638;
t442 = (-t482 * t610 + (-t461 * t606 + t449) * t568) * t638;
t441 = (-t481 * t612 + (t458 * t607 + t446) * t566) * t639;
t440 = (-t480 * t612 + (-t459 * t607 + t447) * t566) * t639;
t439 = (-t469 * t485 + (-t450 * t629 - t462 * t614) * t570) * t494;
t438 = (-t469 * t484 + (-t451 * t629 + t463 * t614) * t570) * t494;
t437 = (-t468 * t483 + (-t448 * t629 - t460 * t615) * t568) * t493;
t436 = (-t468 * t482 + (-t449 * t629 + t461 * t615) * t568) * t493;
t435 = (-t467 * t481 + (-t446 * t629 - t458 * t616) * t566) * t492;
t434 = (-t467 * t480 + (-t447 * t629 + t459 * t616) * t566) * t492;
t433 = (-t454 * t485 + (t450 * t596 - t462 * t617) * t570) * t494;
t432 = (-t454 * t484 + (t451 * t596 + t463 * t617) * t570) * t494;
t431 = (-t453 * t483 + (t448 * t597 - t460 * t618) * t568) * t493;
t430 = (-t453 * t482 + (t449 * t597 + t461 * t618) * t568) * t493;
t429 = (-t452 * t481 + (t446 * t598 - t458 * t619) * t566) * t492;
t428 = (-t452 * t480 + (t447 * t598 + t459 * t619) * t566) * t492;
t1 = [m(4) + (-t433 * t485 + (-t439 * t624 + t445 * t450) * t570) * t494 + (-t431 * t483 + (-t437 * t626 + t443 * t448) * t568) * t493 + (-t429 * t481 + (-t435 * t628 + t441 * t446) * t566) * t492, (-t433 * t484 + (t439 * t623 + t445 * t451) * t570) * t494 + (-t431 * t482 + (t437 * t625 + t443 * t449) * t568) * t493 + (-t429 * t480 + (t435 * t627 + t441 * t447) * t566) * t492, (t439 * t551 + t445 * t497) * t570 + (t437 * t549 + t443 * t496) * t568 + (t435 * t547 + t441 * t495) * t566; (-t432 * t485 + (-t438 * t624 + t444 * t450) * t570) * t494 + (-t430 * t483 + (-t436 * t626 + t442 * t448) * t568) * t493 + (-t428 * t481 + (-t434 * t628 + t440 * t446) * t566) * t492, m(4) + (-t432 * t484 + (t438 * t623 + t444 * t451) * t570) * t494 + (-t430 * t482 + (t436 * t625 + t442 * t449) * t568) * t493 + (-t428 * t480 + (t434 * t627 + t440 * t447) * t566) * t492, (t438 * t551 + t444 * t497) * t570 + (t436 * t549 + t442 * t496) * t568 + (t434 * t547 + t440 * t495) * t566; (-t457 * t485 + (t450 * t475 - t462 * t620) * t570) * t494 + (-t456 * t483 + (t448 * t474 - t460 * t621) * t568) * t493 + (-t455 * t481 + (t446 * t473 - t458 * t622) * t566) * t492, (-t457 * t484 + (t451 * t475 + t463 * t620) * t570) * t494 + (-t456 * t482 + (t449 * t474 + t461 * t621) * t568) * t493 + (-t455 * t480 + (t447 * t473 + t459 * t622) * t566) * t492, m(4) + (t466 * t551 + t475 * t497) * t570 + (t465 * t549 + t474 * t496) * t568 + (t464 * t547 + t473 * t495) * t566;];
MX  = t1;
