% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G3A0
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:26:58
% EndTime: 2020-08-06 19:27:00
% DurationCPUTime: 2.27s
% Computational Cost: add. (6411->313), mult. (7290->534), div. (540->6), fcn. (4050->18), ass. (0->199)
t658 = 2 * rSges(3,3);
t661 = (t658 + qJ(3,3)) * qJ(3,3);
t660 = (t658 + qJ(3,2)) * qJ(3,2);
t659 = (t658 + qJ(3,1)) * qJ(3,1);
t565 = (pkin(2) + pkin(3));
t657 = -2 * t565;
t656 = m(2) * rSges(2,1);
t563 = m(2) * rSges(2,2);
t562 = rSges(2,3) + pkin(5);
t561 = pkin(2) + rSges(3,1);
t548 = sin(qJ(2,3));
t516 = t548 * qJ(3,3);
t554 = cos(qJ(2,3));
t617 = t554 * t565;
t584 = t516 + pkin(1) + t617;
t502 = 0.1e1 / t584;
t655 = m(3) * t502;
t550 = sin(qJ(2,2));
t517 = t550 * qJ(3,2);
t556 = cos(qJ(2,2));
t615 = t556 * t565;
t583 = t517 + pkin(1) + t615;
t503 = 0.1e1 / t583;
t654 = m(3) * t503;
t552 = sin(qJ(2,1));
t518 = t552 * qJ(3,1);
t558 = cos(qJ(2,1));
t613 = t558 * t565;
t582 = t518 + pkin(1) + t613;
t504 = 0.1e1 / t582;
t653 = m(3) * t504;
t560 = pkin(5) + rSges(3,2);
t652 = m(3) * t560;
t651 = m(3) * t561;
t650 = pkin(1) * t548;
t649 = pkin(1) * t550;
t648 = pkin(1) * t552;
t547 = legFrame(1,2);
t521 = sin(t547);
t647 = qJ(3,1) * t521;
t524 = cos(t547);
t646 = qJ(3,1) * t524;
t546 = legFrame(2,2);
t520 = sin(t546);
t645 = qJ(3,2) * t520;
t523 = cos(t546);
t644 = qJ(3,2) * t523;
t545 = legFrame(3,2);
t519 = sin(t545);
t643 = qJ(3,3) * t519;
t522 = cos(t545);
t642 = qJ(3,3) * t522;
t549 = sin(qJ(1,3));
t555 = cos(qJ(1,3));
t564 = pkin(5) - pkin(6);
t611 = t564 * t555;
t641 = (t584 * t549 - t611) * t554;
t551 = sin(qJ(1,2));
t557 = cos(qJ(1,2));
t610 = t564 * t557;
t640 = (t583 * t551 - t610) * t556;
t553 = sin(qJ(1,1));
t559 = cos(qJ(1,1));
t609 = t564 * t559;
t639 = (t582 * t553 - t609) * t558;
t638 = t519 * t549;
t637 = t519 * t565;
t636 = t520 * t551;
t635 = t520 * t565;
t634 = t521 * t553;
t633 = t521 * t565;
t632 = t522 * t549;
t631 = t522 * t565;
t630 = t523 * t551;
t629 = t523 * t565;
t628 = t524 * t553;
t627 = t524 * t565;
t626 = (qJ(3,3) + t565) * (-qJ(3,3) + t565);
t625 = (qJ(3,2) + t565) * (-qJ(3,2) + t565);
t624 = (qJ(3,1) + t565) * (-qJ(3,1) + t565);
t623 = t548 * t560;
t622 = t548 * t565;
t621 = t550 * t560;
t620 = t550 * t565;
t619 = t552 * t560;
t618 = t552 * t565;
t616 = t555 * t565;
t614 = t557 * t565;
t612 = t559 * t565;
t608 = -Icges(2,1) - Icges(3,1);
t607 = pkin(1) * t555 + t549 * t564;
t606 = pkin(1) * t557 + t551 * t564;
t605 = pkin(1) * t559 + t553 * t564;
t604 = rSges(3,3) - t561;
t603 = rSges(3,3) + t561;
t602 = qJ(3,1) * t657;
t601 = qJ(3,2) * t657;
t600 = qJ(3,3) * t657;
t599 = m(3) * t623;
t598 = m(3) * t621;
t597 = m(3) * t619;
t596 = t559 * t518;
t595 = t557 * t517;
t594 = t555 * t516;
t593 = t555 * t626;
t592 = t557 * t625;
t591 = t559 * t624;
t590 = t549 * t623;
t589 = t551 * t621;
t588 = t553 * t619;
t573 = rSges(2,2) ^ 2;
t575 = rSges(2,1) ^ 2;
t587 = Icges(3,2) + Icges(2,3) + (t573 + t575) * m(2);
t572 = rSges(3,3) ^ 2;
t577 = pkin(1) ^ 2;
t586 = t560 ^ 2 + t572 + t577;
t585 = pkin(2) ^ 2 + t572 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t581 = Icges(1,3) + (t562 ^ 2 + t573 + t577) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t608;
t580 = (-t573 + t575) * m(2) + Icges(2,2) + Icges(3,3) + t608;
t579 = -rSges(2,1) * t563 + Icges(2,4) - Icges(3,5);
t578 = -t562 * t563 + Icges(2,6) - Icges(3,6);
t571 = 1 / qJ(3,1);
t569 = 1 / qJ(3,2);
t567 = 1 / qJ(3,3);
t544 = -qJ(3,1) - rSges(3,3);
t543 = -qJ(3,2) - rSges(3,3);
t542 = -qJ(3,3) - rSges(3,3);
t540 = t558 ^ 2;
t539 = t556 ^ 2;
t538 = t554 ^ 2;
t509 = qJ(3,1) + t648;
t508 = qJ(3,2) + t649;
t507 = qJ(3,3) + t650;
t505 = 0.2e1 * pkin(1) * (t651 + t656);
t501 = pkin(1) * qJ(3,1) - t552 * t624;
t500 = pkin(1) * qJ(3,2) - t550 * t625;
t499 = pkin(1) * qJ(3,3) - t548 * t626;
t498 = t560 * t651 + t562 * t656 - Icges(3,4) - Icges(2,5);
t497 = t596 + t605;
t496 = 0.2e1 * t596 + t605;
t495 = t595 + t606;
t494 = 0.2e1 * t595 + t606;
t493 = t594 + t607;
t492 = 0.2e1 * t594 + t607;
t491 = qJ(3,1) * t559 + t605 * t552;
t490 = qJ(3,2) * t557 + t606 * t550;
t489 = qJ(3,3) * t555 + t607 * t548;
t488 = (t585 + t659) * m(3) + t587;
t487 = (t585 + t660) * m(3) + t587;
t486 = (t585 + t661) * m(3) + t587;
t482 = (-t544 * t652 + t578) * t558 - t552 * t498;
t481 = (-t543 * t652 + t578) * t556 - t550 * t498;
t480 = (-t542 * t652 + t578) * t554 - t548 * t498;
t479 = -t553 * t540 * t624 - ((0.2e1 * t518 + pkin(1)) * t553 - t609) * t613 - qJ(3,1) * (t509 * t553 - t552 * t609);
t478 = -t551 * t539 * t625 - ((0.2e1 * t517 + pkin(1)) * t551 - t610) * t615 - qJ(3,2) * (t508 * t551 - t550 * t610);
t477 = -t549 * t538 * t626 - ((0.2e1 * t516 + pkin(1)) * t549 - t611) * t617 - qJ(3,3) * (t507 * t549 - t548 * t611);
t476 = (t524 * t612 - t647) * t540 + (t497 * t524 + t521 * t618) * t558 + t521 * t509;
t475 = (-t521 * t612 - t646) * t540 + (-t497 * t521 + t524 * t618) * t558 + t524 * t509;
t474 = (t523 * t614 - t645) * t539 + (t495 * t523 + t520 * t620) * t556 + t520 * t508;
t473 = (-t520 * t614 - t644) * t539 + (-t495 * t520 + t523 * t620) * t556 + t523 * t508;
t472 = (t522 * t616 - t643) * t538 + (t493 * t522 + t519 * t622) * t554 + t519 * t507;
t471 = (-t519 * t616 - t642) * t538 + (-t493 * t519 + t522 * t622) * t554 + t522 * t507;
t470 = (-(qJ(3,1) + t603) * (qJ(3,1) + t604) * m(3) + t580) * t540 + (0.2e1 * (-t544 * t651 + t579) * t552 + t505) * t558 - 0.2e1 * (m(3) * t544 + t563) * t648 + (t586 + t659) * m(3) + t581;
t469 = (-(qJ(3,2) + t603) * (qJ(3,2) + t604) * m(3) + t580) * t539 + (0.2e1 * (-t543 * t651 + t579) * t550 + t505) * t556 - 0.2e1 * (m(3) * t543 + t563) * t649 + (t586 + t660) * m(3) + t581;
t468 = (-(qJ(3,3) + t603) * (qJ(3,3) + t604) * m(3) + t580) * t538 + (0.2e1 * (-t542 * t651 + t579) * t548 + t505) * t554 - 0.2e1 * (m(3) * t542 + t563) * t650 + (t586 + t661) * m(3) + t581;
t467 = (t521 * t602 + t524 * t591) * t540 + (t496 * t627 - t501 * t521) * t558 + t491 * t646 + t509 * t633;
t466 = (-t521 * t591 + t524 * t602) * t540 + (-t496 * t633 - t501 * t524) * t558 - t491 * t647 + t509 * t627;
t465 = (t520 * t601 + t523 * t592) * t539 + (t494 * t629 - t500 * t520) * t556 + t490 * t644 + t508 * t635;
t464 = (-t520 * t592 + t523 * t601) * t539 + (-t494 * t635 - t500 * t523) * t556 - t490 * t645 + t508 * t629;
t463 = (t519 * t600 + t522 * t593) * t538 + (t492 * t631 - t499 * t519) * t554 + t489 * t642 + t507 * t637;
t462 = (-t519 * t593 + t522 * t600) * t538 + (-t492 * t637 - t499 * t522) * t554 - t489 * t643 + t507 * t631;
t461 = (-t559 * t619 + (t561 * t639 + t479) * t571) * t653;
t460 = (-t557 * t621 + (t561 * t640 + t478) * t569) * t654;
t459 = (-t555 * t623 + (t561 * t641 + t477) * t567) * t655;
t458 = (-t482 * t559 + (-t479 * t651 - t488 * t639) * t571) * t504;
t457 = (-t481 * t557 + (-t478 * t651 - t487 * t640) * t569) * t503;
t456 = (-t480 * t555 + (-t477 * t651 - t486 * t641) * t567) * t502;
t455 = (-t524 * t588 + (-t476 * t561 + t467) * t571) * t653;
t454 = (-t523 * t589 + (-t474 * t561 + t465) * t569) * t654;
t453 = (-t522 * t590 + (-t472 * t561 + t463) * t567) * t655;
t452 = (t521 * t588 + (-t475 * t561 + t466) * t571) * t653;
t451 = (t520 * t589 + (-t473 * t561 + t464) * t569) * t654;
t450 = (t519 * t590 + (-t471 * t561 + t462) * t567) * t655;
t449 = (-t470 * t559 + (t479 * t597 - t482 * t639) * t571) * t504;
t448 = (-t469 * t557 + (t478 * t598 - t481 * t640) * t569) * t503;
t447 = (-t468 * t555 + (t477 * t599 - t480 * t641) * t567) * t502;
t446 = (-t482 * t628 + (-t467 * t651 + t476 * t488) * t571) * t504;
t445 = (-t481 * t630 + (-t465 * t651 + t474 * t487) * t569) * t503;
t444 = (-t480 * t632 + (-t463 * t651 + t472 * t486) * t567) * t502;
t443 = (t482 * t634 + (-t466 * t651 + t475 * t488) * t571) * t504;
t442 = (t481 * t636 + (-t464 * t651 + t473 * t487) * t569) * t503;
t441 = (t480 * t638 + (-t462 * t651 + t471 * t486) * t567) * t502;
t440 = (-t470 * t628 + (t467 * t597 + t476 * t482) * t571) * t504;
t439 = (-t469 * t630 + (t465 * t598 + t474 * t481) * t569) * t503;
t438 = (-t468 * t632 + (t463 * t599 + t472 * t480) * t567) * t502;
t437 = (t470 * t634 + (t466 * t597 + t475 * t482) * t571) * t504;
t436 = (t469 * t636 + (t464 * t598 + t473 * t481) * t569) * t503;
t435 = (t468 * t638 + (t462 * t599 + t471 * t480) * t567) * t502;
t1 = [m(4) + (-t440 * t628 + (t446 * t476 + t455 * t467) * t571) * t504 + (-t439 * t630 + (t445 * t474 + t454 * t465) * t569) * t503 + (-t438 * t632 + (t444 * t472 + t453 * t463) * t567) * t502, (t440 * t634 + (t446 * t475 + t455 * t466) * t571) * t504 + (t439 * t636 + (t445 * t473 + t454 * t464) * t569) * t503 + (t438 * t638 + (t444 * t471 + t453 * t462) * t567) * t502, (-t440 * t559 + (-t446 * t639 + t455 * t479) * t571) * t504 + (-t439 * t557 + (-t445 * t640 + t454 * t478) * t569) * t503 + (-t438 * t555 + (-t444 * t641 + t453 * t477) * t567) * t502; (-t437 * t628 + (t443 * t476 + t452 * t467) * t571) * t504 + (-t436 * t630 + (t442 * t474 + t451 * t465) * t569) * t503 + (-t435 * t632 + (t441 * t472 + t450 * t463) * t567) * t502, m(4) + (t437 * t634 + (t443 * t475 + t452 * t466) * t571) * t504 + (t436 * t636 + (t442 * t473 + t451 * t464) * t569) * t503 + (t435 * t638 + (t441 * t471 + t450 * t462) * t567) * t502, (-t437 * t559 + (-t443 * t639 + t452 * t479) * t571) * t504 + (-t436 * t557 + (-t442 * t640 + t451 * t478) * t569) * t503 + (-t435 * t555 + (-t441 * t641 + t450 * t477) * t567) * t502; (-t449 * t628 + (t458 * t476 + t461 * t467) * t571) * t504 + (-t448 * t630 + (t457 * t474 + t460 * t465) * t569) * t503 + (-t447 * t632 + (t456 * t472 + t459 * t463) * t567) * t502, (t449 * t634 + (t458 * t475 + t461 * t466) * t571) * t504 + (t448 * t636 + (t457 * t473 + t460 * t464) * t569) * t503 + (t447 * t638 + (t456 * t471 + t459 * t462) * t567) * t502, m(4) + (-t449 * t559 + (-t458 * t639 + t461 * t479) * t571) * t504 + (-t448 * t557 + (-t457 * t640 + t460 * t478) * t569) * t503 + (-t447 * t555 + (-t456 * t641 + t459 * t477) * t567) * t502;];
MX  = t1;
