% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:02:30
% EndTime: 2020-08-06 21:02:31
% DurationCPUTime: 1.36s
% Computational Cost: add. (3573->250), mult. (4938->384), div. (249->9), fcn. (2013->44), ass. (0->179)
t686 = rSges(2,3) + pkin(5);
t685 = 2 * pkin(1);
t684 = m(2) / 0.2e1;
t683 = Icges(2,2) / 0.2e1;
t682 = Icges(3,2) / 0.2e1;
t641 = pkin(2) ^ 2;
t681 = t641 / 0.2e1;
t680 = -2 * pkin(1);
t678 = m(2) * rSges(2,1);
t677 = m(2) * rSges(2,2);
t676 = m(3) * rSges(3,1);
t602 = (qJ(2,3) + pkin(7));
t564 = 2 * t602;
t559 = sin(t564);
t633 = 2 * qJ(2,3);
t601 = t633 + pkin(7);
t573 = sin(t601);
t574 = sin(t602);
t608 = sin(t633);
t616 = sin(qJ(2,3));
t640 = pkin(3) ^ 2;
t651 = 2 * pkin(2) * pkin(3);
t521 = t573 * t651 + t640 * t559 + t641 * t608 + (pkin(2) * t616 + pkin(3) * t574) * t685;
t675 = t521 / 0.2e1;
t604 = qJ(2,2) + pkin(7);
t565 = 2 * t604;
t560 = sin(t565);
t634 = 2 * qJ(2,2);
t603 = t634 + pkin(7);
t575 = sin(t603);
t576 = sin(t604);
t609 = sin(t634);
t618 = sin(qJ(2,2));
t522 = t575 * t651 + t640 * t560 + t641 * t609 + (pkin(2) * t618 + pkin(3) * t576) * t685;
t674 = t522 / 0.2e1;
t606 = qJ(2,1) + pkin(7);
t566 = 2 * t606;
t561 = sin(t566);
t635 = 2 * qJ(2,1);
t605 = t635 + pkin(7);
t577 = sin(t605);
t578 = sin(t606);
t610 = sin(t635);
t620 = sin(qJ(2,1));
t523 = t577 * t651 + t640 * t561 + t641 * t610 + (pkin(2) * t620 + pkin(3) * t578) * t685;
t673 = t523 / 0.2e1;
t637 = rSges(2,2) ^ 2;
t639 = rSges(2,1) ^ 2;
t672 = m(3) * t681 + (-t637 + t639) * t684 + t683 - Icges(2,1) / 0.2e1;
t636 = rSges(3,2) ^ 2;
t638 = rSges(3,1) ^ 2;
t671 = (-t636 + t638) * m(3) / 0.2e1 - Icges(3,1) / 0.2e1 + t682;
t622 = cos(qJ(2,3));
t591 = t622 * pkin(2);
t568 = t591 + pkin(1);
t579 = cos(t602);
t525 = -t579 * rSges(3,1) + t574 * rSges(3,2) - t568;
t670 = m(3) * t525;
t624 = cos(qJ(2,2));
t592 = t624 * pkin(2);
t569 = t592 + pkin(1);
t580 = cos(t604);
t526 = -t580 * rSges(3,1) + t576 * rSges(3,2) - t569;
t669 = m(3) * t526;
t626 = cos(qJ(2,1));
t593 = t626 * pkin(2);
t570 = t593 + pkin(1);
t581 = cos(t606);
t527 = -t581 * rSges(3,1) + t578 * rSges(3,2) - t570;
t668 = m(3) * t527;
t655 = pkin(5) + qJ(3,3);
t597 = -pkin(6) - t655;
t588 = 0.1e1 / t597;
t667 = m(3) * t588;
t656 = pkin(5) + qJ(3,2);
t598 = -pkin(6) - t656;
t589 = 0.1e1 / t598;
t666 = m(3) * t589;
t657 = pkin(5) + qJ(3,1);
t599 = -pkin(6) - t657;
t590 = 0.1e1 / t599;
t665 = m(3) * t590;
t594 = rSges(3,3) + t655;
t664 = m(3) * t594;
t595 = rSges(3,3) + t656;
t663 = m(3) * t595;
t596 = rSges(3,3) + t657;
t662 = m(3) * t596;
t661 = pkin(3) * t579;
t660 = pkin(3) * t580;
t659 = pkin(3) * t581;
t611 = sin(pkin(7));
t658 = pkin(3) * t611;
t612 = cos(pkin(7));
t563 = t612 * pkin(3) + pkin(2);
t654 = 0.1e1 / (t563 * t622 - t616 * t658) * (t616 * t563 + t622 * t658);
t653 = 0.1e1 / (t563 * t624 - t618 * t658) * (t618 * t563 + t624 * t658);
t652 = 0.1e1 / (t563 * t626 - t620 * t658) * (t620 * t563 + t626 * t658);
t650 = t637 + t639;
t648 = t588 * t654;
t647 = t589 * t653;
t646 = t590 * t652;
t645 = t686 * t677 - Icges(2,6);
t644 = -t686 * t678 + Icges(2,5);
t632 = 2 * pkin(1) ^ 2;
t643 = t632 / 0.2e1 + t636 / 0.2e1 + t638 / 0.2e1 + t681;
t562 = t612 * pkin(2) * t676;
t642 = Icges(1,3) + ((2 * pkin(5) ^ 2) + t632 + ((4 * pkin(5) + 2 * rSges(2,3)) * rSges(2,3)) + t650) * t684 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t562 + t682 + t683 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t627 = cos(qJ(1,1));
t625 = cos(qJ(1,2));
t623 = cos(qJ(1,3));
t621 = sin(qJ(1,1));
t619 = sin(qJ(1,2));
t617 = sin(qJ(1,3));
t615 = legFrame(1,3);
t614 = legFrame(2,3);
t613 = legFrame(3,3);
t587 = cos(t615);
t586 = cos(t614);
t585 = cos(t613);
t584 = sin(t615);
t583 = sin(t614);
t582 = sin(t613);
t572 = -rSges(2,1) * t677 + Icges(2,4);
t571 = -rSges(3,2) * t676 + Icges(3,4);
t567 = pkin(2) * m(3) + t678;
t558 = rSges(3,1) * t662 - Icges(3,5);
t557 = rSges(3,1) * t663 - Icges(3,5);
t556 = rSges(3,1) * t664 - Icges(3,5);
t555 = rSges(3,2) * t662 - Icges(3,6);
t554 = rSges(3,2) * t663 - Icges(3,6);
t553 = rSges(3,2) * t664 - Icges(3,6);
t549 = 0.1e1 / (t593 + t659);
t548 = 0.1e1 / (t592 + t660);
t547 = 0.1e1 / (t591 + t661);
t545 = t584 * t627 + t587 * t621;
t544 = t583 * t625 + t586 * t619;
t543 = t582 * t623 + t585 * t617;
t542 = -t584 * t621 + t587 * t627;
t541 = -t583 * t619 + t586 * t625;
t540 = -t582 * t617 + t585 * t623;
t539 = t570 * t627 - t621 * t599;
t538 = t569 * t625 - t619 * t598;
t537 = t568 * t623 - t617 * t597;
t536 = t621 * t570 + t627 * t599;
t535 = t619 * t569 + t625 * t598;
t534 = t617 * t568 + t623 * t597;
t524 = 0.2e1 * t562 + t650 * m(2) + Icges(2,3) + Icges(3,3) + (-0.2e1 * t611 * rSges(3,2) * pkin(2) + t636 + t638 + t641) * m(3);
t520 = t536 * t587 + t539 * t584 + t545 * t659;
t519 = t535 * t586 + t538 * t583 + t544 * t660;
t518 = t534 * t585 + t537 * t582 + t543 * t661;
t517 = -t536 * t584 + t539 * t587 + t542 * t659;
t516 = -t535 * t583 + t538 * t586 + t541 * t660;
t515 = -t534 * t582 + t537 * t585 + t540 * t661;
t514 = (-pkin(2) * t662 + t555 * t611 - t558 * t612 + t644) * t620 - (t555 * t612 + t558 * t611 + t645) * t626;
t513 = (-pkin(2) * t663 + t554 * t611 - t557 * t612 + t644) * t618 - (t554 * t612 + t557 * t611 + t645) * t624;
t512 = (-pkin(2) * t664 + t553 * t611 - t556 * t612 + t644) * t616 - (t553 * t612 + t556 * t611 + t645) * t622;
t511 = (t527 * t652 + t549 * t673) * t665;
t510 = (t526 * t653 + t548 * t674) * t666;
t509 = (t525 * t654 + t547 * t675) * t667;
t508 = (t527 * t545 + t520) * t665;
t507 = (t526 * t544 + t519) * t666;
t506 = (t525 * t543 + t518) * t667;
t505 = (t527 * t542 + t517) * t665;
t504 = (t526 * t541 + t516) * t666;
t503 = (t525 * t540 + t515) * t667;
t502 = cos(t566) * t671 + cos(t635) * t672 + (t567 * t626 - t620 * t677) * t685 + t571 * t561 + t572 * t610 + (t596 ^ 2 + (pkin(2) * cos(t605) + t581 * t685) * rSges(3,1) + (t578 * t680 + (-t577 - t611) * pkin(2)) * rSges(3,2) + t643) * m(3) + t642;
t501 = cos(t565) * t671 + cos(t634) * t672 + (t567 * t624 - t618 * t677) * t685 + t571 * t560 + t572 * t609 + (t595 ^ 2 + (pkin(2) * cos(t603) + t580 * t685) * rSges(3,1) + (t576 * t680 + (-t575 - t611) * pkin(2)) * rSges(3,2) + t643) * m(3) + t642;
t500 = cos(t564) * t671 + cos(t633) * t672 + (t567 * t622 - t616 * t677) * t685 + t571 * t559 + t572 * t608 + (t594 ^ 2 + (pkin(2) * cos(t601) + t579 * t685) * rSges(3,1) + (t574 * t680 + (-t573 - t611) * pkin(2)) * rSges(3,2) + t643) * m(3) + t642;
t499 = (t502 * t545 + t520 * t668) * t590;
t498 = (t501 * t544 + t519 * t669) * t589;
t497 = (t500 * t543 + t518 * t670) * t588;
t496 = (t502 * t542 + t517 * t668) * t590;
t495 = (t501 * t541 + t516 * t669) * t589;
t494 = (t500 * t540 + t515 * t670) * t588;
t493 = -t502 * t646 + (t514 - t523 * t527 * t665 / 0.2e1) * t549;
t492 = -t501 * t647 + (t513 - t522 * t526 * t666 / 0.2e1) * t548;
t491 = -t500 * t648 + (t512 - t521 * t525 * t667 / 0.2e1) * t547;
t1 = [m(4) - (-t496 * t542 - t505 * t517) * t590 - (-t495 * t541 - t504 * t516) * t589 - (-t494 * t540 - t503 * t515) * t588, -(-t496 * t545 - t505 * t520) * t590 - (-t495 * t544 - t504 * t519) * t589 - (-t494 * t543 - t503 * t518) * t588, -(-t496 * t652 + (-t505 * t673 + t542 * t514) * t549) * t590 - (-t495 * t653 + (-t504 * t674 + t541 * t513) * t548) * t589 - (-t494 * t654 + (-t503 * t675 + t540 * t512) * t547) * t588; -(-t499 * t542 - t508 * t517) * t590 - (-t498 * t541 - t507 * t516) * t589 - (-t497 * t540 - t506 * t515) * t588, m(4) - (-t499 * t545 - t508 * t520) * t590 - (-t498 * t544 - t507 * t519) * t589 - (-t497 * t543 - t506 * t518) * t588, -(-t499 * t652 + (-t508 * t673 + t545 * t514) * t549) * t590 - (-t498 * t653 + (-t507 * t674 + t544 * t513) * t548) * t589 - (-t497 * t654 + (-t506 * t675 + t543 * t512) * t547) * t588; -(t493 * t542 - t511 * t517) * t590 - (t492 * t541 - t510 * t516) * t589 - (t491 * t540 - t509 * t515) * t588, -(t493 * t545 - t511 * t520) * t590 - (t492 * t544 - t510 * t519) * t589 - (t491 * t543 - t509 * t518) * t588, -t491 * t648 - t492 * t647 - t493 * t646 + m(4) + (t549 * t524 - (-t511 * t673 + t514 * t652) * t590) * t549 + (t548 * t524 - (-t510 * t674 + t513 * t653) * t589) * t548 + (t547 * t524 - (-t509 * t675 + t512 * t654) * t588) * t547;];
MX  = t1;
