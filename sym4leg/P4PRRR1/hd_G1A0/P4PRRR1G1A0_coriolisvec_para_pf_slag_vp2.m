% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:12:32
% EndTime: 2020-03-02 20:12:39
% DurationCPUTime: 6.97s
% Computational Cost: add. (21829->251), mult. (14176->507), div. (3248->6), fcn. (15388->34), ass. (0->276)
t709 = 0.1e1 / pkin(2);
t678 = pkin(7) + qJ(2,1);
t664 = qJ(3,1) + t678;
t645 = sin(t664);
t648 = cos(t664);
t658 = sin(t678);
t661 = cos(t678);
t816 = 0.1e1 / (t645 * t661 - t648 * t658);
t827 = t816 * t709;
t677 = pkin(7) + qJ(2,2);
t663 = qJ(3,2) + t677;
t644 = sin(t663);
t647 = cos(t663);
t657 = sin(t677);
t660 = cos(t677);
t817 = 0.1e1 / (t644 * t660 - t647 * t657);
t826 = t817 * t709;
t676 = pkin(7) + qJ(2,3);
t662 = qJ(3,3) + t676;
t643 = sin(t662);
t646 = cos(t662);
t656 = sin(t676);
t659 = cos(t676);
t818 = 0.1e1 / (t643 * t659 - t646 * t656);
t825 = t818 * t709;
t675 = pkin(7) + qJ(2,4);
t655 = qJ(3,4) + t675;
t641 = sin(t655);
t642 = cos(t655);
t653 = sin(t675);
t654 = cos(t675);
t819 = 0.1e1 / (t641 * t654 - t642 * t653);
t824 = t819 * t709;
t695 = xP(4);
t673 = sin(t695);
t674 = cos(t695);
t701 = koppelP(1,2);
t705 = koppelP(1,1);
t632 = -t673 * t701 + t674 * t705;
t692 = xDP(4);
t693 = xDP(2);
t616 = t632 * t692 + t693;
t628 = t673 * t705 + t674 * t701;
t694 = xDP(1);
t620 = -t628 * t692 + t694;
t683 = legFrame(1,3);
t668 = sin(t683);
t672 = cos(t683);
t612 = -t645 * t668 + t648 * t672;
t745 = t612 * t827;
t611 = t645 * t672 + t648 * t668;
t746 = t611 * t827;
t566 = t616 * t746 + t620 * t745;
t707 = 0.1e1 / pkin(3);
t769 = t707 * t709;
t739 = t816 * t769;
t592 = -pkin(2) * (t658 * t668 - t661 * t672) + t612 * pkin(3);
t798 = t592 * t620;
t589 = pkin(2) * (t658 * t672 + t661 * t668) + t611 * pkin(3);
t801 = t589 * t616;
t533 = -(t798 / 0.2e1 + t801 / 0.2e1) * t739 + t566;
t710 = (t798 + t801) * t739;
t536 = -t710 + t566;
t706 = pkin(3) ^ 2;
t708 = pkin(2) ^ 2;
t718 = t645 * t658 + t648 * t661;
t714 = t718 * pkin(2);
t763 = 0.2e1 * pkin(2) * pkin(3);
t810 = t536 * t710;
t520 = ((-t533 * t718 * t763 - t536 * t706 - t708 * t566) * t707 * t566 + (pkin(3) + t714) * t810) * t827;
t524 = (-pkin(3) * t810 + (pkin(3) * t536 + t566 * t714) * t566) * t827;
t688 = sin(qJ(3,1));
t814 = mrSges(3,2) * pkin(2);
t652 = t688 * t814;
t691 = cos(qJ(3,1));
t815 = mrSges(3,1) * pkin(2);
t764 = t691 * t815;
t636 = Ifges(3,3) - t652 + t764;
t516 = Ifges(3,3) * t520 + t524 * t636;
t771 = t816 * t707;
t640 = mrSges(3,1) * t688 + mrSges(3,2) * t691;
t806 = t566 ^ 2 * t640;
t823 = -t516 * t739 - t771 * t806;
t700 = koppelP(2,2);
t704 = koppelP(2,1);
t631 = -t673 * t700 + t674 * t704;
t615 = t631 * t692 + t693;
t627 = t673 * t704 + t674 * t700;
t619 = -t627 * t692 + t694;
t682 = legFrame(2,3);
t667 = sin(t682);
t671 = cos(t682);
t610 = -t644 * t667 + t647 * t671;
t747 = t610 * t826;
t609 = t644 * t671 + t647 * t667;
t748 = t609 * t826;
t565 = t615 * t748 + t619 * t747;
t741 = t817 * t769;
t591 = -pkin(2) * (t657 * t667 - t660 * t671) + t610 * pkin(3);
t799 = t591 * t619;
t588 = pkin(2) * (t657 * t671 + t660 * t667) + t609 * pkin(3);
t802 = t588 * t615;
t532 = -(t799 / 0.2e1 + t802 / 0.2e1) * t741 + t565;
t711 = (t799 + t802) * t741;
t535 = -t711 + t565;
t719 = t644 * t657 + t647 * t660;
t715 = t719 * pkin(2);
t811 = t535 * t711;
t519 = ((-t532 * t719 * t763 - t535 * t706 - t708 * t565) * t707 * t565 + (pkin(3) + t715) * t811) * t826;
t523 = (-pkin(3) * t811 + (pkin(3) * t535 + t565 * t715) * t565) * t826;
t687 = sin(qJ(3,2));
t651 = t687 * t814;
t690 = cos(qJ(3,2));
t765 = t690 * t815;
t635 = Ifges(3,3) - t651 + t765;
t515 = Ifges(3,3) * t519 + t523 * t635;
t773 = t817 * t707;
t639 = mrSges(3,1) * t687 + mrSges(3,2) * t690;
t807 = t565 ^ 2 * t639;
t822 = -t515 * t741 - t773 * t807;
t699 = koppelP(3,2);
t703 = koppelP(3,1);
t630 = -t673 * t699 + t674 * t703;
t614 = t630 * t692 + t693;
t626 = t673 * t703 + t674 * t699;
t618 = -t626 * t692 + t694;
t681 = legFrame(3,3);
t666 = sin(t681);
t670 = cos(t681);
t608 = -t643 * t666 + t646 * t670;
t749 = t608 * t825;
t607 = t643 * t670 + t646 * t666;
t750 = t607 * t825;
t564 = t614 * t750 + t618 * t749;
t743 = t818 * t769;
t590 = -pkin(2) * (t656 * t666 - t659 * t670) + t608 * pkin(3);
t800 = t590 * t618;
t587 = pkin(2) * (t656 * t670 + t659 * t666) + t607 * pkin(3);
t803 = t587 * t614;
t531 = -(t800 / 0.2e1 + t803 / 0.2e1) * t743 + t564;
t712 = (t800 + t803) * t743;
t534 = -t712 + t564;
t720 = t643 * t656 + t646 * t659;
t716 = t720 * pkin(2);
t812 = t534 * t712;
t518 = ((-t531 * t720 * t763 - t534 * t706 - t708 * t564) * t707 * t564 + (pkin(3) + t716) * t812) * t825;
t522 = (-pkin(3) * t812 + (pkin(3) * t534 + t564 * t716) * t564) * t825;
t686 = sin(qJ(3,3));
t650 = t686 * t814;
t689 = cos(qJ(3,3));
t766 = t689 * t815;
t634 = Ifges(3,3) - t650 + t766;
t514 = Ifges(3,3) * t518 + t522 * t634;
t775 = t818 * t707;
t638 = mrSges(3,1) * t686 + mrSges(3,2) * t689;
t808 = t564 ^ 2 * t638;
t821 = -t514 * t743 - t775 * t808;
t698 = koppelP(4,2);
t702 = koppelP(4,1);
t629 = -t673 * t698 + t674 * t702;
t613 = t629 * t692 + t693;
t625 = t673 * t702 + t674 * t698;
t617 = -t625 * t692 + t694;
t680 = legFrame(4,3);
t665 = sin(t680);
t669 = cos(t680);
t606 = -t641 * t665 + t642 * t669;
t753 = t606 * t824;
t605 = t641 * t669 + t642 * t665;
t754 = t605 * t824;
t558 = t613 * t754 + t617 * t753;
t751 = t819 * t769;
t586 = -pkin(2) * (t653 * t665 - t669 * t654) + t606 * pkin(3);
t804 = t586 * t617;
t585 = pkin(2) * (t653 * t669 + t654 * t665) + t605 * pkin(3);
t805 = t585 * t613;
t526 = -(t804 / 0.2e1 + t805 / 0.2e1) * t751 + t558;
t713 = (t804 + t805) * t751;
t527 = -t713 + t558;
t721 = t641 * t653 + t642 * t654;
t717 = t721 * pkin(2);
t813 = t527 * t713;
t517 = ((-t526 * t721 * t763 - t527 * t706 - t708 * t558) * t707 * t558 + (pkin(3) + t717) * t813) * t824;
t521 = (-pkin(3) * t813 + (pkin(3) * t527 + t558 * t717) * t558) * t824;
t684 = sin(qJ(3,4));
t649 = t684 * t814;
t685 = cos(qJ(3,4));
t767 = t685 * t815;
t633 = Ifges(3,3) - t649 + t767;
t510 = Ifges(3,3) * t517 + t521 * t633;
t792 = t819 * t707;
t637 = mrSges(3,1) * t684 + mrSges(3,2) * t685;
t809 = t558 ^ 2 * t637;
t820 = -t510 * t751 - t792 * t809;
t679 = t692 ^ 2;
t797 = t819 * t605;
t796 = t819 * t606;
t726 = m(3) * t708 + Ifges(2,3) + Ifges(3,3);
t621 = -0.2e1 * t649 + t726 + 0.2e1 * t767;
t795 = t819 * t621;
t794 = t819 * t633;
t790 = t818 * t607;
t789 = t818 * t608;
t622 = -0.2e1 * t650 + t726 + 0.2e1 * t766;
t788 = t818 * t622;
t787 = t818 * t634;
t785 = t817 * t609;
t784 = t817 * t610;
t623 = -0.2e1 * t651 + t726 + 0.2e1 * t765;
t783 = t817 * t623;
t782 = t817 * t635;
t780 = t816 * t611;
t779 = t816 * t612;
t624 = -0.2e1 * t652 + t726 + 0.2e1 * t764;
t778 = t816 * t624;
t777 = t816 * t636;
t768 = t709 * t679;
t762 = t585 * t792;
t761 = t586 * t792;
t760 = t587 * t775;
t759 = t588 * t773;
t758 = t589 * t771;
t757 = t590 * t775;
t756 = t591 * t773;
t755 = t592 * t771;
t752 = t633 * t792;
t744 = t634 * t775;
t742 = t635 * t773;
t740 = t636 * t771;
t738 = 0.2e1 * t526 * t713 * t637;
t737 = 0.2e1 * t531 * t712 * t638;
t736 = 0.2e1 * t532 * t711 * t639;
t735 = 0.2e1 * t533 * t710 * t640;
t725 = t819 * t738;
t724 = t818 * t737;
t723 = t817 * t736;
t722 = t816 * t735;
t697 = mrSges(4,1);
t696 = mrSges(4,2);
t576 = (t611 * t632 - t612 * t628) * t827;
t575 = (t609 * t631 - t610 * t627) * t826;
t574 = (t607 * t630 - t608 * t626) * t825;
t573 = (t605 * t629 - t606 * t625) * t824;
t572 = (-Ifges(3,3) * t755 + t612 * t777) * t709;
t571 = (-Ifges(3,3) * t758 + t611 * t777) * t709;
t570 = (-Ifges(3,3) * t756 + t610 * t782) * t709;
t569 = (-Ifges(3,3) * t759 + t609 * t782) * t709;
t568 = (-Ifges(3,3) * t757 + t608 * t787) * t709;
t567 = (-Ifges(3,3) * t760 + t607 * t787) * t709;
t560 = (-Ifges(3,3) * t761 + t606 * t794) * t709;
t559 = (-Ifges(3,3) * t762 + t605 * t794) * t709;
t556 = (-t592 * t740 + t612 * t778) * t709;
t555 = (-t589 * t740 + t611 * t778) * t709;
t554 = (-t591 * t742 + t610 * t783) * t709;
t553 = (-t588 * t742 + t609 * t783) * t709;
t552 = (-t590 * t744 + t608 * t788) * t709;
t551 = (-t587 * t744 + t607 * t788) * t709;
t550 = (-t586 * t752 + t606 * t795) * t709;
t549 = (-t585 * t752 + t605 * t795) * t709;
t548 = (t589 * t632 - t592 * t628) * t739;
t547 = (t588 * t631 - t591 * t627) * t741;
t546 = (t587 * t630 - t590 * t626) * t743;
t545 = (t585 * t629 - t586 * t625) * t751;
t540 = -Ifges(3,3) * t548 + t576 * t636;
t539 = -Ifges(3,3) * t547 + t575 * t635;
t538 = -Ifges(3,3) * t546 + t574 * t634;
t537 = -Ifges(3,3) * t545 + t573 * t633;
t530 = -t548 * t636 + t576 * t624;
t529 = -t547 * t635 + t575 * t623;
t528 = -t546 * t634 + t574 * t622;
t525 = -t545 * t633 + t573 * t621;
t513 = t520 * t636 + t524 * t624;
t512 = t519 * t635 + t523 * t623;
t511 = t518 * t634 + t522 * t622;
t509 = t517 * t633 + t521 * t621;
t1 = [t513 * t745 + t612 * t722 + t512 * t747 + t610 * t723 + t511 * t749 + t608 * t724 + t509 * t753 + t606 * t725 + t679 * (t673 * t696 - t674 * t697) + t823 * t592 + t822 * t591 + t821 * t590 + t820 * t586 + (-(t556 * t779 - t572 * t755) * t632 - (t556 * t780 - t572 * t758) * t628 - (t554 * t784 - t570 * t756) * t631 - (t554 * t785 - t570 * t759) * t627 - (t552 * t789 - t568 * t757) * t630 - (t552 * t790 - t568 * t760) * t626 - (t550 * t796 - t560 * t761) * t629 - (t550 * t797 - t560 * t762) * t625) * t768; t513 * t746 + t611 * t722 + t512 * t748 + t609 * t723 + t511 * t750 + t607 * t724 + t509 * t754 + t605 * t725 - t679 * (t673 * t697 + t674 * t696) + t823 * t589 + t822 * t588 + t821 * t587 + t820 * t585 + (-(t555 * t779 - t571 * t755) * t632 - (t555 * t780 - t571 * t758) * t628 - (t553 * t784 - t569 * t756) * t631 - (t553 * t785 - t569 * t759) * t627 - (t551 * t789 - t567 * t757) * t630 - (t551 * t790 - t567 * t760) * t626 - (t549 * t796 - t559 * t761) * t629 - (t549 * t797 - t559 * t762) * t625) * t768; 0; t573 * t509 - t545 * t510 + t574 * t511 + t575 * t512 + t576 * t513 - t546 * t514 - t547 * t515 - t548 * t516 + (-(t530 * t779 - t540 * t755) * t632 - (t530 * t780 - t540 * t758) * t628 - (t529 * t784 - t539 * t756) * t631 - (t529 * t785 - t539 * t759) * t627 - (t528 * t789 - t538 * t757) * t630 - (t528 * t790 - t538 * t760) * t626 - (t525 * t796 - t537 * t761) * t629 - (t525 * t797 - t537 * t762) * t625) * t768 + (-t545 * t809 - t546 * t808 - t547 * t807 - t548 * t806 + t573 * t738 + t574 * t737 + t575 * t736 + t576 * t735) * pkin(2);];
taucX  = t1;
