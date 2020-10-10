% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:31
% EndTime: 2020-08-06 19:20:33
% DurationCPUTime: 2.29s
% Computational Cost: add. (6411->313), mult. (7290->533), div. (540->6), fcn. (4050->18), ass. (0->201)
t658 = 2 * rSges(3,3);
t661 = (t658 + qJ(3,3)) * qJ(3,3);
t660 = (t658 + qJ(3,2)) * qJ(3,2);
t659 = (t658 + qJ(3,1)) * qJ(3,1);
t563 = (pkin(2) + pkin(3));
t657 = -2 * t563;
t656 = m(2) * rSges(2,1);
t561 = m(2) * rSges(2,2);
t560 = rSges(2,3) + pkin(5);
t559 = pkin(2) + rSges(3,1);
t546 = sin(qJ(2,3));
t514 = t546 * qJ(3,3);
t552 = cos(qJ(2,3));
t609 = t563 * t552;
t582 = t514 + pkin(1) + t609;
t500 = 0.1e1 / t582;
t655 = m(3) * t500;
t548 = sin(qJ(2,2));
t515 = t548 * qJ(3,2);
t554 = cos(qJ(2,2));
t608 = t563 * t554;
t581 = t515 + pkin(1) + t608;
t501 = 0.1e1 / t581;
t654 = m(3) * t501;
t550 = sin(qJ(2,1));
t516 = t550 * qJ(3,1);
t556 = cos(qJ(2,1));
t607 = t563 * t556;
t580 = t516 + pkin(1) + t607;
t502 = 0.1e1 / t580;
t653 = m(3) * t502;
t652 = pkin(1) * t546;
t651 = pkin(1) * t548;
t650 = pkin(1) * t550;
t540 = -qJ(3,3) - rSges(3,3);
t649 = t540 * m(3);
t541 = -qJ(3,2) - rSges(3,3);
t648 = t541 * m(3);
t542 = -qJ(3,1) - rSges(3,3);
t647 = t542 * m(3);
t646 = t559 * m(3);
t547 = sin(qJ(1,3));
t562 = pkin(5) - pkin(6);
t508 = t547 * t562;
t553 = cos(qJ(1,3));
t645 = (t582 * t553 + t508) * t552;
t549 = sin(qJ(1,2));
t509 = t549 * t562;
t555 = cos(qJ(1,2));
t644 = (t581 * t555 + t509) * t554;
t551 = sin(qJ(1,1));
t510 = t551 * t562;
t557 = cos(qJ(1,1));
t643 = (t580 * t557 + t510) * t556;
t586 = pkin(1) * t547 - t562 * t553;
t622 = t547 * qJ(3,3);
t593 = t546 * t622;
t642 = (t586 + 0.2e1 * t593) * t563;
t585 = pkin(1) * t549 - t562 * t555;
t618 = t549 * qJ(3,2);
t591 = t548 * t618;
t641 = (t585 + 0.2e1 * t591) * t563;
t584 = pkin(1) * t551 - t562 * t557;
t614 = t551 * qJ(3,1);
t589 = t550 * t614;
t640 = (t584 + 0.2e1 * t589) * t563;
t543 = legFrame(3,2);
t517 = sin(t543);
t639 = t517 * qJ(3,3);
t638 = t517 * t553;
t544 = legFrame(2,2);
t518 = sin(t544);
t637 = t518 * qJ(3,2);
t636 = t518 * t555;
t545 = legFrame(1,2);
t519 = sin(t545);
t635 = t519 * qJ(3,1);
t634 = t519 * t557;
t520 = cos(t543);
t633 = t520 * qJ(3,3);
t632 = t520 * t553;
t521 = cos(t544);
t631 = t521 * qJ(3,2);
t630 = t521 * t555;
t522 = cos(t545);
t629 = t522 * qJ(3,1);
t628 = t522 * t557;
t627 = (qJ(3,3) + t563) * (-qJ(3,3) + t563);
t626 = (qJ(3,2) + t563) * (-qJ(3,2) + t563);
t625 = (qJ(3,1) + t563) * (-qJ(3,1) + t563);
t558 = pkin(5) + rSges(3,2);
t624 = t546 * t558;
t623 = t546 * t563;
t621 = t547 * t563;
t620 = t548 * t558;
t619 = t548 * t563;
t617 = t549 * t563;
t616 = t550 * t558;
t615 = t550 * t563;
t613 = t551 * t563;
t505 = qJ(3,3) + t652;
t612 = t563 * t505;
t506 = qJ(3,2) + t651;
t611 = t563 * t506;
t507 = qJ(3,1) + t650;
t610 = t563 * t507;
t606 = -Icges(2,1) - Icges(3,1);
t605 = rSges(3,3) - t559;
t604 = rSges(3,3) + t559;
t603 = qJ(3,1) * t657;
t602 = qJ(3,2) * t657;
t601 = qJ(3,3) * t657;
t600 = m(3) * t624;
t599 = m(3) * t620;
t598 = m(3) * t616;
t597 = t547 * t627;
t596 = t549 * t626;
t595 = t551 * t625;
t594 = t553 * t624;
t592 = t555 * t620;
t590 = t557 * t616;
t571 = rSges(2,2) ^ 2;
t573 = rSges(2,1) ^ 2;
t588 = Icges(3,2) + Icges(2,3) + (t571 + t573) * m(2);
t570 = rSges(3,3) ^ 2;
t575 = pkin(1) ^ 2;
t587 = t558 ^ 2 + t570 + t575;
t583 = pkin(2) ^ 2 + t570 + (2 * pkin(2) + rSges(3,1)) * rSges(3,1);
t579 = Icges(1,3) + (t560 ^ 2 + t571 + t575) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - t606;
t578 = (-t571 + t573) * m(2) + Icges(2,2) + Icges(3,3) + t606;
t577 = -rSges(2,1) * t561 + Icges(2,4) - Icges(3,5);
t576 = -t560 * t561 + Icges(2,6) - Icges(3,6);
t569 = 1 / qJ(3,1);
t567 = 1 / qJ(3,2);
t565 = 1 / qJ(3,3);
t538 = t556 ^ 2;
t537 = t554 ^ 2;
t536 = t552 ^ 2;
t503 = 0.2e1 * (t646 + t656) * pkin(1);
t499 = pkin(1) * qJ(3,1) - t550 * t625;
t498 = pkin(1) * qJ(3,2) - t548 * t626;
t497 = pkin(1) * qJ(3,3) - t546 * t627;
t496 = t558 * t646 + t560 * t656 - Icges(3,4) - Icges(2,5);
t495 = t584 + t589;
t493 = t585 + t591;
t491 = t586 + t593;
t489 = t584 * t550 + t614;
t488 = t585 * t548 + t618;
t487 = t586 * t546 + t622;
t486 = (t583 + t659) * m(3) + t588;
t485 = (t583 + t660) * m(3) + t588;
t484 = (t583 + t661) * m(3) + t588;
t480 = (-t558 * t647 + t576) * t556 - t550 * t496;
t479 = (-t558 * t648 + t576) * t554 - t548 * t496;
t478 = (-t558 * t649 + t576) * t552 - t546 * t496;
t477 = t557 * t538 * t625 + ((0.2e1 * t516 + pkin(1)) * t557 + t510) * t607 + qJ(3,1) * (t507 * t557 + t550 * t510);
t476 = t555 * t537 * t626 + ((0.2e1 * t515 + pkin(1)) * t555 + t509) * t608 + qJ(3,2) * (t506 * t555 + t548 * t509);
t475 = t553 * t536 * t627 + ((0.2e1 * t514 + pkin(1)) * t553 + t508) * t609 + qJ(3,3) * (t505 * t553 + t546 * t508);
t474 = (t522 * t613 - t635) * t538 + (t495 * t522 + t519 * t615) * t556 + t519 * t507;
t473 = (-t519 * t613 - t629) * t538 + (-t495 * t519 + t522 * t615) * t556 + t522 * t507;
t472 = (t521 * t617 - t637) * t537 + (t493 * t521 + t518 * t619) * t554 + t518 * t506;
t471 = (-t518 * t617 - t631) * t537 + (-t493 * t518 + t521 * t619) * t554 + t521 * t506;
t470 = (t520 * t621 - t639) * t536 + (t491 * t520 + t517 * t623) * t552 + t517 * t505;
t469 = (-t517 * t621 - t633) * t536 + (-t491 * t517 + t520 * t623) * t552 + t520 * t505;
t468 = (-(qJ(3,1) + t604) * (qJ(3,1) + t605) * m(3) + t578) * t538 + (0.2e1 * (-t542 * t646 + t577) * t550 + t503) * t556 - 0.2e1 * (t561 + t647) * t650 + (t587 + t659) * m(3) + t579;
t467 = (-(qJ(3,2) + t604) * (qJ(3,2) + t605) * m(3) + t578) * t537 + (0.2e1 * (-t541 * t646 + t577) * t548 + t503) * t554 - 0.2e1 * (t561 + t648) * t651 + (t587 + t660) * m(3) + t579;
t466 = (-(qJ(3,3) + t604) * (qJ(3,3) + t605) * m(3) + t578) * t536 + (0.2e1 * (-t540 * t646 + t577) * t546 + t503) * t552 - 0.2e1 * (t561 + t649) * t652 + (t587 + t661) * m(3) + t579;
t465 = (t519 * t603 + t522 * t595) * t538 + (-t519 * t499 + t522 * t640) * t556 + t489 * t629 + t519 * t610;
t464 = (-t519 * t595 + t522 * t603) * t538 + (-t522 * t499 - t519 * t640) * t556 - t489 * t635 + t522 * t610;
t463 = (t518 * t602 + t521 * t596) * t537 + (-t518 * t498 + t521 * t641) * t554 + t488 * t631 + t518 * t611;
t462 = (-t518 * t596 + t521 * t602) * t537 + (-t521 * t498 - t518 * t641) * t554 - t488 * t637 + t521 * t611;
t461 = (t517 * t601 + t520 * t597) * t536 + (-t517 * t497 + t520 * t642) * t552 + t487 * t633 + t517 * t612;
t460 = (-t517 * t597 + t520 * t601) * t536 + (-t520 * t497 - t517 * t642) * t552 - t487 * t639 + t520 * t612;
t459 = (-t551 * t616 + (-t559 * t643 + t477) * t569) * t653;
t458 = (-t549 * t620 + (-t559 * t644 + t476) * t567) * t654;
t457 = (-t547 * t624 + (-t559 * t645 + t475) * t565) * t655;
t456 = (-t480 * t551 + (-t477 * t646 + t486 * t643) * t569) * t502;
t455 = (-t479 * t549 + (-t476 * t646 + t485 * t644) * t567) * t501;
t454 = (-t478 * t547 + (-t475 * t646 + t484 * t645) * t565) * t500;
t453 = (t522 * t590 + (-t474 * t559 + t465) * t569) * t653;
t452 = (t521 * t592 + (-t472 * t559 + t463) * t567) * t654;
t451 = (t520 * t594 + (-t470 * t559 + t461) * t565) * t655;
t450 = (-t519 * t590 + (-t473 * t559 + t464) * t569) * t653;
t449 = (-t518 * t592 + (-t471 * t559 + t462) * t567) * t654;
t448 = (-t517 * t594 + (-t469 * t559 + t460) * t565) * t655;
t447 = (-t468 * t551 + (t477 * t598 + t480 * t643) * t569) * t502;
t446 = (-t467 * t549 + (t476 * t599 + t479 * t644) * t567) * t501;
t445 = (-t466 * t547 + (t475 * t600 + t478 * t645) * t565) * t500;
t444 = (t480 * t628 + (-t465 * t646 + t474 * t486) * t569) * t502;
t443 = (t479 * t630 + (-t463 * t646 + t472 * t485) * t567) * t501;
t442 = (t478 * t632 + (-t461 * t646 + t470 * t484) * t565) * t500;
t441 = (-t480 * t634 + (-t464 * t646 + t473 * t486) * t569) * t502;
t440 = (-t479 * t636 + (-t462 * t646 + t471 * t485) * t567) * t501;
t439 = (-t478 * t638 + (-t460 * t646 + t469 * t484) * t565) * t500;
t438 = (t468 * t628 + (t465 * t598 + t474 * t480) * t569) * t502;
t437 = (t467 * t630 + (t463 * t599 + t472 * t479) * t567) * t501;
t436 = (t466 * t632 + (t461 * t600 + t470 * t478) * t565) * t500;
t435 = (-t468 * t634 + (t464 * t598 + t473 * t480) * t569) * t502;
t434 = (-t467 * t636 + (t462 * t599 + t471 * t479) * t567) * t501;
t433 = (-t466 * t638 + (t460 * t600 + t469 * t478) * t565) * t500;
t1 = [m(4) + (t438 * t628 + (t444 * t474 + t453 * t465) * t569) * t502 + (t437 * t630 + (t443 * t472 + t452 * t463) * t567) * t501 + (t436 * t632 + (t442 * t470 + t451 * t461) * t565) * t500, (-t438 * t634 + (t444 * t473 + t453 * t464) * t569) * t502 + (-t437 * t636 + (t443 * t471 + t452 * t462) * t567) * t501 + (-t436 * t638 + (t442 * t469 + t451 * t460) * t565) * t500, (-t438 * t551 + (t444 * t643 + t453 * t477) * t569) * t502 + (-t437 * t549 + (t443 * t644 + t452 * t476) * t567) * t501 + (-t436 * t547 + (t442 * t645 + t451 * t475) * t565) * t500; (t435 * t628 + (t441 * t474 + t450 * t465) * t569) * t502 + (t434 * t630 + (t440 * t472 + t449 * t463) * t567) * t501 + (t433 * t632 + (t439 * t470 + t448 * t461) * t565) * t500, m(4) + (-t435 * t634 + (t441 * t473 + t450 * t464) * t569) * t502 + (-t434 * t636 + (t440 * t471 + t449 * t462) * t567) * t501 + (-t433 * t638 + (t439 * t469 + t448 * t460) * t565) * t500, (-t435 * t551 + (t441 * t643 + t450 * t477) * t569) * t502 + (-t434 * t549 + (t440 * t644 + t449 * t476) * t567) * t501 + (-t433 * t547 + (t439 * t645 + t448 * t475) * t565) * t500; (t447 * t628 + (t456 * t474 + t459 * t465) * t569) * t502 + (t446 * t630 + (t455 * t472 + t458 * t463) * t567) * t501 + (t445 * t632 + (t454 * t470 + t457 * t461) * t565) * t500, (-t447 * t634 + (t456 * t473 + t459 * t464) * t569) * t502 + (-t446 * t636 + (t455 * t471 + t458 * t462) * t567) * t501 + (-t445 * t638 + (t454 * t469 + t457 * t460) * t565) * t500, m(4) + (-t447 * t551 + (t456 * t643 + t459 * t477) * t569) * t502 + (-t446 * t549 + (t455 * t644 + t458 * t476) * t567) * t501 + (-t445 * t547 + (t454 * t645 + t457 * t475) * t565) * t500;];
MX  = t1;
