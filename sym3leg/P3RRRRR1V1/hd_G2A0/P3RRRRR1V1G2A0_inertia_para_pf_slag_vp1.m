% Calculate inertia matrix for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:31
% EndTime: 2020-08-07 03:33:33
% DurationCPUTime: 1.44s
% Computational Cost: add. (3846->265), mult. (6768->440), div. (648->8), fcn. (3474->54), ass. (0->224)
t764 = 2 * pkin(1);
t763 = m(2) / 0.2e1;
t762 = m(3) / 0.2e1;
t761 = Icges(2,2) / 0.2e1;
t760 = Icges(3,2) / 0.2e1;
t759 = -2 * pkin(1);
t757 = pkin(2) * m(3);
t756 = m(2) * rSges(2,2);
t755 = m(2) * rSges(2,3);
t754 = m(3) * rSges(3,2);
t753 = m(3) * rSges(3,3);
t632 = qJ(2,3) + qJ(3,3);
t752 = cos(qJ(1,3) - t632) / 0.2e1 + cos(qJ(1,3) + t632) / 0.2e1;
t634 = qJ(2,2) + qJ(3,2);
t751 = cos(qJ(1,2) - t634) / 0.2e1 + cos(qJ(1,2) + t634) / 0.2e1;
t636 = qJ(2,1) + qJ(3,1);
t750 = cos(qJ(1,1) - t636) / 0.2e1 + cos(qJ(1,1) + t636) / 0.2e1;
t667 = rSges(2,2) ^ 2;
t669 = rSges(2,1) ^ 2;
t671 = pkin(2) ^ 2;
t749 = t671 * t762 + (-t667 + t669) * t763 + t761 - Icges(2,1) / 0.2e1;
t666 = rSges(3,2) ^ 2;
t668 = rSges(3,1) ^ 2;
t748 = (-t666 + t668) * t762 - Icges(3,1) / 0.2e1 + t760;
t637 = legFrame(3,2);
t625 = cos(t637);
t649 = cos(qJ(3,3));
t650 = cos(qJ(2,3));
t640 = sin(qJ(3,3));
t641 = sin(qJ(2,3));
t710 = t640 * t641;
t678 = -t649 * t650 + t710;
t709 = t640 * t650;
t679 = -t641 * t649 - t709;
t622 = sin(t637);
t642 = sin(qJ(1,3));
t719 = t622 * t642;
t549 = t679 * t625 + t678 * t719;
t628 = 0.1e1 / t640;
t747 = t549 * t628;
t716 = t625 * t642;
t550 = t679 * t622 - t678 * t716;
t746 = t550 * t628;
t638 = legFrame(2,2);
t626 = cos(t638);
t652 = cos(qJ(3,2));
t653 = cos(qJ(2,2));
t643 = sin(qJ(3,2));
t644 = sin(qJ(2,2));
t708 = t643 * t644;
t676 = -t652 * t653 + t708;
t707 = t643 * t653;
t677 = -t644 * t652 - t707;
t623 = sin(t638);
t645 = sin(qJ(1,2));
t718 = t623 * t645;
t551 = t677 * t626 + t676 * t718;
t629 = 0.1e1 / t643;
t745 = t551 * t629;
t715 = t626 * t645;
t552 = t677 * t623 - t676 * t715;
t744 = t552 * t629;
t639 = legFrame(1,2);
t627 = cos(t639);
t655 = cos(qJ(3,1));
t656 = cos(qJ(2,1));
t646 = sin(qJ(3,1));
t647 = sin(qJ(2,1));
t706 = t646 * t647;
t674 = -t655 * t656 + t706;
t705 = t646 * t656;
t675 = -t647 * t655 - t705;
t624 = sin(t639);
t648 = sin(qJ(1,1));
t717 = t624 * t648;
t553 = t675 * t627 + t674 * t717;
t630 = 0.1e1 / t646;
t743 = t553 * t630;
t714 = t627 * t648;
t554 = t675 * t624 - t674 * t714;
t742 = t554 * t630;
t603 = pkin(3) * t649 + pkin(2);
t573 = pkin(3) * t709 + t641 * t603;
t576 = -pkin(3) * t710 + t603 * t650;
t555 = t622 * t573 - t576 * t716;
t741 = t555 * t628;
t604 = pkin(3) * t652 + pkin(2);
t574 = pkin(3) * t707 + t644 * t604;
t577 = -pkin(3) * t708 + t604 * t653;
t556 = t623 * t574 - t577 * t715;
t740 = t556 * t629;
t605 = pkin(3) * t655 + pkin(2);
t575 = pkin(3) * t705 + t647 * t605;
t578 = -pkin(3) * t706 + t605 * t656;
t557 = t624 * t575 - t578 * t714;
t739 = t557 * t630;
t558 = t573 * t625 + t576 * t719;
t738 = t558 * t628;
t559 = t574 * t626 + t577 * t718;
t737 = t559 * t629;
t560 = t575 * t627 + t578 * t717;
t736 = t560 * t630;
t611 = rSges(3,2) * t753 - Icges(3,6);
t612 = rSges(3,1) * t753 - Icges(3,5);
t616 = sin(t632);
t619 = cos(t632);
t567 = t619 * t611 + t612 * t616;
t670 = 0.1e1 / pkin(3);
t735 = t567 * t670;
t617 = sin(t634);
t620 = cos(t634);
t568 = t620 * t611 + t612 * t617;
t734 = t568 * t670;
t618 = sin(t636);
t621 = cos(t636);
t569 = t621 * t611 + t612 * t618;
t733 = t569 * t670;
t704 = t666 + t668;
t597 = t704 * m(3) + Icges(3,3);
t701 = rSges(3,1) * t757;
t599 = t649 * t701;
t700 = pkin(2) * t754;
t683 = t640 * t700;
t570 = t599 - t683 + t597;
t732 = t570 * t670;
t600 = t652 * t701;
t682 = t643 * t700;
t571 = t600 - t682 + t597;
t731 = t571 * t670;
t601 = t655 * t701;
t681 = t646 * t700;
t572 = t601 - t681 + t597;
t730 = t572 * t670;
t651 = cos(qJ(1,3));
t729 = t576 * t651;
t654 = cos(qJ(1,2));
t728 = t577 * t654;
t657 = cos(qJ(1,1));
t727 = t578 * t657;
t579 = 0.1e1 / (t650 * pkin(2) + pkin(3) * t619 + pkin(1));
t726 = t579 * t642;
t725 = t579 * t651;
t580 = 0.1e1 / (t653 * pkin(2) + pkin(3) * t620 + pkin(1));
t724 = t580 * t645;
t723 = t580 * t654;
t581 = 0.1e1 / (t656 * pkin(2) + pkin(3) * t621 + pkin(1));
t722 = t581 * t648;
t721 = t581 * t657;
t720 = t597 * t670;
t672 = 0.1e1 / pkin(2);
t713 = t628 * t672;
t712 = t629 * t672;
t711 = t630 * t672;
t703 = t667 + t669;
t699 = t628 * t729;
t698 = t670 * t729;
t697 = t629 * t728;
t696 = t670 * t728;
t695 = t630 * t727;
t694 = t670 * t727;
t693 = t622 * t725;
t692 = t625 * t725;
t691 = t623 * t723;
t690 = t626 * t723;
t689 = t624 * t721;
t688 = t627 * t721;
t687 = t671 + t704;
t686 = t628 * t752;
t685 = t629 * t751;
t684 = t630 * t750;
t680 = t703 * m(2) + t687 * m(3) + Icges(2,3) + Icges(3,3);
t662 = 2 * pkin(1) ^ 2;
t673 = Icges(1,3) + (0.2e1 * rSges(3,3) ^ 2 + t662 + t687) * t762 + (0.2e1 * rSges(2,3) ^ 2 + t662 + t703) * t763 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t760 + t761 + Icges(3,1) / 0.2e1 + Icges(2,1) / 0.2e1;
t665 = 0.2e1 * qJ(2,1);
t664 = 0.2e1 * qJ(2,2);
t663 = 0.2e1 * qJ(2,3);
t635 = t665 + qJ(3,1);
t633 = qJ(3,2) + t664;
t631 = t663 + qJ(3,3);
t615 = -rSges(2,1) * t756 + Icges(2,4);
t614 = rSges(2,2) * t755 - Icges(2,6);
t613 = -rSges(3,1) * t754 + Icges(3,4);
t610 = 0.2e1 * t636;
t609 = 0.2e1 * t634;
t608 = 0.2e1 * t632;
t602 = m(2) * rSges(2,1) + t757;
t594 = rSges(2,1) * t755 + pkin(2) * t753 - Icges(2,5);
t566 = 0.2e1 * t601 + t680 - 0.2e1 * t681;
t565 = 0.2e1 * t600 + t680 - 0.2e1 * t682;
t564 = 0.2e1 * t599 + t680 - 0.2e1 * t683;
t563 = t594 * t647 + t656 * t614 + t569;
t562 = t594 * t644 + t653 * t614 + t568;
t561 = t594 * t641 + t650 * t614 + t567;
t548 = -t569 * t722 + (t572 * t750 - t597 * t694) * t711;
t547 = -t568 * t724 + (t571 * t751 - t597 * t696) * t712;
t546 = -t567 * t726 + (t570 * t752 - t597 * t698) * t713;
t545 = -t563 * t722 + (t566 * t750 - t572 * t694) * t711;
t544 = -t562 * t724 + (t565 * t751 - t571 * t696) * t712;
t543 = -t561 * t726 + (t564 * t752 - t570 * t698) * t713;
t542 = cos(t610) * t748 + cos(t665) * t749 + t601 + (t602 * t656 - t647 * t756) * t764 + t613 * sin(t610) + t615 * sin(t665) + ((cos(t635) * pkin(2) + t621 * t764) * rSges(3,1) + (t618 * t759 + (-sin(t635) - t646) * pkin(2)) * rSges(3,2)) * m(3) + t673;
t541 = cos(t609) * t748 + cos(t664) * t749 + t600 + (t602 * t653 - t644 * t756) * t764 + t613 * sin(t609) + t615 * sin(t664) + ((cos(t633) * pkin(2) + t620 * t764) * rSges(3,1) + (t617 * t759 + (-sin(t633) - t643) * pkin(2)) * rSges(3,2)) * m(3) + t673;
t540 = cos(t608) * t748 + cos(t663) * t749 + t599 + (t602 * t650 - t641 * t756) * t764 + t613 * sin(t608) + t615 * sin(t663) + ((cos(t631) * pkin(2) + t619 * t764) * rSges(3,1) + (t616 * t759 + (-sin(t631) - t640) * pkin(2)) * rSges(3,2)) * m(3) + t673;
t539 = t569 * t688 + (t554 * t572 + t557 * t720) * t711;
t538 = t568 * t690 + (t552 * t571 + t556 * t720) * t712;
t537 = t567 * t692 + (t550 * t570 + t555 * t720) * t713;
t536 = -t569 * t689 + (t553 * t572 + t560 * t720) * t711;
t535 = -t568 * t691 + (t551 * t571 + t559 * t720) * t712;
t534 = -t567 * t693 + (t549 * t570 + t558 * t720) * t713;
t533 = t563 * t688 + (t554 * t566 + t557 * t730) * t711;
t532 = t562 * t690 + (t552 * t565 + t556 * t731) * t712;
t531 = t561 * t692 + (t550 * t564 + t555 * t732) * t713;
t530 = -t563 * t689 + (t553 * t566 + t560 * t730) * t711;
t529 = -t562 * t691 + (t551 * t565 + t559 * t731) * t712;
t528 = -t561 * t693 + (t549 * t564 + t558 * t732) * t713;
t527 = -t542 * t722 + (t563 * t750 - t569 * t694) * t711;
t526 = -t541 * t724 + (t562 * t751 - t568 * t696) * t712;
t525 = -t540 * t726 + (t561 * t752 - t567 * t698) * t713;
t524 = t542 * t688 + (t554 * t563 + t557 * t733) * t711;
t523 = t541 * t690 + (t552 * t562 + t556 * t734) * t712;
t522 = t540 * t692 + (t550 * t561 + t555 * t735) * t713;
t521 = -t542 * t689 + (t553 * t563 + t560 * t733) * t711;
t520 = -t541 * t691 + (t551 * t562 + t559 * t734) * t712;
t519 = -t540 * t693 + (t549 * t561 + t558 * t735) * t713;
t1 = [t522 * t692 + t523 * t690 + t524 * t688 + m(4) + (t531 * t746 + t532 * t744 + t533 * t742 + (t537 * t741 + t538 * t740 + t539 * t739) * t670) * t672, -t522 * t693 - t523 * t691 - t524 * t689 + (t531 * t747 + t532 * t745 + t533 * t743 + (t537 * t738 + t538 * t737 + t539 * t736) * t670) * t672, -t522 * t726 - t523 * t724 - t524 * t722 + (t533 * t684 + t532 * t685 + t531 * t686 + (-t537 * t699 - t538 * t697 - t539 * t695) * t670) * t672; t519 * t692 + t520 * t690 + t521 * t688 + (t528 * t746 + t529 * t744 + t530 * t742 + (t534 * t741 + t535 * t740 + t536 * t739) * t670) * t672, -t519 * t693 - t520 * t691 - t521 * t689 + m(4) + (t528 * t747 + t529 * t745 + t530 * t743 + (t534 * t738 + t535 * t737 + t536 * t736) * t670) * t672, -t519 * t726 - t520 * t724 - t521 * t722 + (t530 * t684 + t529 * t685 + t528 * t686 + (-t534 * t699 - t535 * t697 - t536 * t695) * t670) * t672; t525 * t692 + t526 * t690 + t527 * t688 + (t543 * t746 + t544 * t744 + t545 * t742 + (t546 * t741 + t547 * t740 + t548 * t739) * t670) * t672, -t525 * t693 - t526 * t691 - t527 * t689 + (t543 * t747 + t544 * t745 + t545 * t743 + (t546 * t738 + t547 * t737 + t548 * t736) * t670) * t672, -t525 * t726 - t526 * t724 - t527 * t722 + m(4) + (t545 * t684 + t544 * t685 + t543 * t686 + (-t546 * t699 - t547 * t697 - t548 * t695) * t670) * t672;];
MX  = t1;
