% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:16
% EndTime: 2020-08-06 19:59:19
% DurationCPUTime: 2.87s
% Computational Cost: add. (13980->270), mult. (23631->467), div. (3051->12), fcn. (20190->35), ass. (0->212)
t630 = pkin(4) + qJ(3,3);
t639 = sin(qJ(1,3));
t645 = cos(qJ(1,3));
t629 = cos(pkin(5));
t729 = pkin(2) * t629;
t591 = pkin(1) + t729;
t644 = cos(qJ(2,3));
t628 = sin(pkin(5));
t638 = sin(qJ(2,3));
t715 = t628 * t638;
t668 = pkin(2) * t715 - t591 * t644;
t558 = -t630 * t645 - t668 * t639;
t730 = pkin(2) * t628;
t568 = t591 * t638 + t644 * t730;
t635 = legFrame(3,2);
t609 = sin(t635);
t612 = cos(t635);
t530 = t558 * t612 + t568 * t609;
t531 = -t558 * t609 + t568 * t612;
t615 = t644 * pkin(1);
t619 = qJ(2,3) + pkin(5);
t756 = t615 + pkin(2) * cos(t619);
t562 = t630 * t639 + t645 * t756;
t622 = 0.1e1 / t630;
t651 = xDP(3);
t652 = xDP(2);
t653 = xDP(1);
t521 = (t530 * t653 + t531 * t652 + t562 * t651) * t622;
t765 = 0.2e1 * t521;
t631 = pkin(4) + qJ(3,2);
t641 = sin(qJ(1,2));
t647 = cos(qJ(1,2));
t646 = cos(qJ(2,2));
t640 = sin(qJ(2,2));
t714 = t628 * t640;
t667 = pkin(2) * t714 - t591 * t646;
t559 = -t631 * t647 - t667 * t641;
t569 = t591 * t640 + t646 * t730;
t636 = legFrame(2,2);
t610 = sin(t636);
t613 = cos(t636);
t532 = t559 * t613 + t569 * t610;
t533 = -t559 * t610 + t569 * t613;
t616 = t646 * pkin(1);
t620 = qJ(2,2) + pkin(5);
t755 = t616 + pkin(2) * cos(t620);
t563 = t631 * t641 + t647 * t755;
t623 = 0.1e1 / t631;
t522 = (t532 * t653 + t533 * t652 + t563 * t651) * t623;
t764 = 0.2e1 * t522;
t632 = pkin(4) + qJ(3,1);
t643 = sin(qJ(1,1));
t649 = cos(qJ(1,1));
t648 = cos(qJ(2,1));
t642 = sin(qJ(2,1));
t713 = t628 * t642;
t666 = pkin(2) * t713 - t591 * t648;
t560 = -t632 * t649 - t666 * t643;
t570 = t591 * t642 + t648 * t730;
t637 = legFrame(1,2);
t611 = sin(t637);
t614 = cos(t637);
t534 = t560 * t614 + t570 * t611;
t535 = -t560 * t611 + t570 * t614;
t617 = t648 * pkin(1);
t621 = qJ(2,1) + pkin(5);
t754 = t617 + pkin(2) * cos(t621);
t564 = t632 * t643 + t649 * t754;
t624 = 0.1e1 / t632;
t523 = (t534 * t653 + t535 * t652 + t564 * t651) * t624;
t763 = 0.2e1 * t523;
t707 = 0.2e1 * pkin(1);
t741 = 0.2e1 * t629;
t762 = t741 / 0.2e1;
t760 = mrSges(3,2) * t629;
t574 = 0.1e1 / t756;
t548 = (t609 * t653 + t612 * t652) * t574;
t545 = t548 ^ 2;
t659 = pkin(2) ^ 2;
t660 = pkin(1) ^ 2;
t708 = -t659 - t660;
t585 = t729 * t707 - t708;
t759 = t545 * t585;
t575 = 0.1e1 / t755;
t549 = (t610 * t653 + t613 * t652) * t575;
t546 = t549 ^ 2;
t758 = t546 * t585;
t576 = 0.1e1 / t754;
t550 = (t611 * t653 + t614 * t652) * t576;
t547 = t550 ^ 2;
t757 = t547 * t585;
t728 = mrSges(3,1) * t628;
t753 = pkin(1) * t728 - Ifges(2,4) + Ifges(3,4);
t601 = mrSges(3,2) * qJ(3,1) - Ifges(3,6);
t604 = mrSges(3,1) * qJ(3,1) - Ifges(3,5);
t739 = m(3) * qJ(3,1);
t608 = mrSges(3,3) + t739;
t752 = -pkin(1) * t608 + t601 * t628 - t604 * t629 + Ifges(2,5);
t600 = mrSges(3,2) * qJ(3,2) - Ifges(3,6);
t603 = mrSges(3,1) * qJ(3,2) - Ifges(3,5);
t738 = m(3) * qJ(3,2);
t607 = mrSges(3,3) + t738;
t751 = -pkin(1) * t607 + t600 * t628 - t603 * t629 + Ifges(2,5);
t599 = mrSges(3,2) * qJ(3,3) - Ifges(3,6);
t602 = mrSges(3,1) * qJ(3,3) - Ifges(3,5);
t737 = m(3) * qJ(3,3);
t606 = mrSges(3,3) + t737;
t750 = -pkin(1) * t606 + t599 * t628 - t602 * t629 + Ifges(2,5);
t633 = Ifges(3,2) - Ifges(3,1);
t654 = m(3) * t660;
t727 = mrSges(3,2) * t628;
t572 = -0.2e1 * pkin(1) * t727 - Ifges(2,1) + Ifges(2,2) - t633 + t654;
t618 = t629 ^ 2;
t712 = t633 * t618;
t726 = Ifges(3,4) * t628;
t740 = pkin(1) * mrSges(3,1);
t551 = (0.4e1 * t726 + 0.2e1 * t740) * t629 + t572 + 0.2e1 * t712;
t749 = -0.2e1 * mrSges(3,3);
t699 = t639 * t730;
t718 = t591 * t639;
t536 = (-t609 * t718 + t612 * t730) * t644 + t638 * (t591 * t612 + t609 * t699);
t539 = (t609 * t730 + t612 * t718) * t644 + (t591 * t609 - t612 * t699) * t638;
t565 = 0.1e1 / t668;
t518 = (t645 * t651 - (t536 * t652 + t539 * t653) * t565) * t622;
t748 = 0.2e1 * t518;
t698 = t641 * t730;
t717 = t591 * t641;
t537 = (-t610 * t717 + t613 * t730) * t646 + t640 * (t591 * t613 + t610 * t698);
t540 = (t610 * t730 + t613 * t717) * t646 + (t591 * t610 - t613 * t698) * t640;
t566 = 0.1e1 / t667;
t519 = (t647 * t651 - (t537 * t652 + t540 * t653) * t566) * t623;
t747 = 0.2e1 * t519;
t697 = t643 * t730;
t716 = t591 * t643;
t538 = (-t611 * t716 + t614 * t730) * t648 + t642 * (t591 * t614 + t611 * t697);
t541 = (t611 * t730 + t614 * t716) * t648 + (t591 * t611 - t614 * t697) * t642;
t567 = 0.1e1 / t666;
t520 = (t649 * t651 - (t538 * t652 + t541 * t653) * t567) * t624;
t746 = 0.2e1 * t520;
t745 = -0.2e1 * t521;
t744 = -0.2e1 * t522;
t743 = -0.2e1 * t523;
t584 = pkin(1) * mrSges(3,2) + t628 * t633;
t725 = t618 * Ifges(3,4);
t561 = t584 * t629 - 0.2e1 * t725 + t753;
t742 = 0.2e1 * t561;
t724 = t518 * t638;
t723 = t519 * t640;
t722 = t520 * t642;
t721 = t584 * t518;
t720 = t584 * t519;
t719 = t584 * t520;
t577 = t638 * pkin(1) + pkin(2) * sin(t619);
t656 = 0.2e1 * qJ(2,3);
t693 = pkin(2) * t707;
t500 = (-t759 + ((-t659 * cos(0.2e1 * t619) - cos(t656) * t660 + (-cos(t656 + pkin(5)) - t629) * t693 + t708) * t518 / 0.2e1 + t756 * t765 + (-t518 * t630 + 0.2e1 * t548 * t577) * t630) * t518) * t622;
t512 = (-0.1e1 / (t615 + (t629 * t644 - t715) * pkin(2)) * t759 + (t668 * t518 + t765) * t518) * t622;
t686 = -t602 * t628 + Ifges(2,6);
t542 = -t750 * t638 + (t599 * t629 - t686) * t644;
t655 = m(3) * pkin(1);
t672 = -mrSges(3,1) * t629 + t727;
t573 = -t655 + t672;
t583 = t728 + t760;
t552 = -t573 * t644 - t583 * t638;
t625 = t644 ^ 2;
t705 = 0.2e1 * t726;
t661 = t629 * t705 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t712;
t662 = 0.4e1 * (t705 + t740) * t629 + 0.2e1 * t572 + 0.4e1 * t712;
t689 = t545 * t574 * t577;
t692 = t518 * t561 * t625;
t696 = t518 * t725;
t711 = (t638 * t644 * t742 - t551 * t625 + (t749 - t737) * qJ(3,3) + t661) * t512 - t542 * t689 + t552 * t500 + t521 * t606 * t748 - t545 * t686 * t638 + (-t662 * t724 * t644 + t753 * t748 + t721 * t741 - 0.4e1 * t692 - 0.4e1 * t696 + (t599 * t638 * t762 + t750 * t644) * t548) * t548;
t578 = t640 * pkin(1) + pkin(2) * sin(t620);
t657 = 0.2e1 * qJ(2,2);
t501 = (-t758 + ((-t659 * cos(0.2e1 * t620) - cos(t657) * t660 + (-cos(pkin(5) + t657) - t629) * t693 + t708) * t519 / 0.2e1 + t755 * t764 + (-t519 * t631 + 0.2e1 * t549 * t578) * t631) * t519) * t623;
t513 = (-0.1e1 / (t616 + (t629 * t646 - t714) * pkin(2)) * t758 + (t667 * t519 + t764) * t519) * t623;
t685 = -t603 * t628 + Ifges(2,6);
t543 = -t751 * t640 + (t600 * t629 - t685) * t646;
t553 = -t573 * t646 - t583 * t640;
t626 = t646 ^ 2;
t688 = t546 * t575 * t578;
t691 = t519 * t561 * t626;
t695 = t519 * t725;
t710 = (t640 * t646 * t742 - t551 * t626 + (t749 - t738) * qJ(3,2) + t661) * t513 - t543 * t688 + t553 * t501 + t522 * t607 * t747 - t546 * t685 * t640 + (-t662 * t723 * t646 + t753 * t747 + t720 * t741 - 0.4e1 * t691 - 0.4e1 * t695 + (t600 * t640 * t762 + t751 * t646) * t549) * t549;
t579 = pkin(1) * t642 + pkin(2) * sin(t621);
t658 = 0.2e1 * qJ(2,1);
t502 = (-t757 + ((-t659 * cos(0.2e1 * t621) - cos(t658) * t660 + (-cos(pkin(5) + t658) - t629) * t693 + t708) * t520 / 0.2e1 + t754 * t763 + (-t520 * t632 + 0.2e1 * t550 * t579) * t632) * t520) * t624;
t514 = (-0.1e1 / (t617 + (t629 * t648 - t713) * pkin(2)) * t757 + (t666 * t520 + t763) * t520) * t624;
t684 = -t604 * t628 + Ifges(2,6);
t544 = -t752 * t642 + (t601 * t629 - t684) * t648;
t554 = -t573 * t648 - t583 * t642;
t627 = t648 ^ 2;
t687 = t547 * t576 * t579;
t690 = t520 * t561 * t627;
t694 = t520 * t725;
t709 = (t642 * t648 * t742 - t551 * t627 + (t749 - t739) * qJ(3,1) + t661) * t514 - t544 * t687 + t554 * t502 + t523 * t608 * t746 - t547 * t684 * t642 + (-t662 * t722 * t648 + t753 * t746 + t719 * t741 - 0.4e1 * t690 - 0.4e1 * t694 + (t601 * t642 * t762 + t752 * t648) * t550) * t550;
t706 = -0.2e1 * t728;
t703 = t638 * t745;
t702 = t640 * t744;
t701 = t642 * t743;
t683 = t711 * t565;
t682 = t710 * t566;
t681 = t709 * t567;
t680 = -m(3) * t500 + t512 * t552 + (-t518 * t606 / 0.2e1 + (-t573 * t638 + t583 * t644) * t548) * t748;
t679 = -m(3) * t501 + t513 * t553 + (-t519 * t607 / 0.2e1 + (-t573 * t640 + t583 * t646) * t549) * t747;
t678 = -m(3) * t502 + t514 * t554 + (-t520 * t608 / 0.2e1 + (-t573 * t642 + t583 * t648) * t550) * t746;
t571 = t672 * t707 - Ifges(2,3) - Ifges(3,3) - t654;
t589 = t655 - t727;
t671 = t574 * ((0.2e1 * t692 + (t521 * t706 + t551 * t724 + t745 * t760) * t644 + 0.2e1 * t696 + (mrSges(3,1) * t703 - t721) * t629 + t589 * t703 - t518 * t753) * t518 + t512 * t542 - t571 * t689);
t670 = t575 * ((0.2e1 * t691 + (t522 * t706 + t551 * t723 + t744 * t760) * t646 + 0.2e1 * t695 + (mrSges(3,1) * t702 - t720) * t629 + t589 * t702 - t519 * t753) * t519 + t513 * t543 - t571 * t688);
t669 = t576 * ((0.2e1 * t690 + (t523 * t706 + t551 * t722 + t743 * t760) * t648 + 0.2e1 * t694 + (mrSges(3,1) * t701 - t719) * t629 + t589 * t701 - t520 * t753) * t520 + t514 * t544 - t571 * t687);
t1 = [t611 * t669 + t610 * t670 + t609 * t671 + (t678 * t534 - t541 * t681) * t624 + (t679 * t532 - t540 * t682) * t623 + (t680 * t530 - t539 * t683) * t622; t614 * t669 + t613 * t670 + t612 * t671 + (t678 * t535 - t538 * t681) * t624 + (t679 * t533 - t537 * t682) * t623 + (t680 * t531 - t536 * t683) * t622; (t678 * t564 + t709 * t649) * t624 + (t679 * t563 + t710 * t647) * t623 + (t680 * t562 + t711 * t645) * t622;];
taucX  = t1;
