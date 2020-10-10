% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x13]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:34
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:34:04
% EndTime: 2020-08-06 19:34:07
% DurationCPUTime: 2.21s
% Computational Cost: add. (7031->239), mult. (11310->555), div. (1548->21), fcn. (8961->18), ass. (0->225)
t850 = 2 * pkin(1);
t709 = pkin(1) ^ 2;
t849 = -2 * t709;
t699 = cos(qJ(2,3));
t669 = 0.1e1 / t699;
t693 = sin(qJ(2,3));
t820 = t669 * t693;
t701 = cos(qJ(2,2));
t674 = 0.1e1 / t701;
t695 = sin(qJ(2,2));
t819 = t674 * t695;
t703 = cos(qJ(2,1));
t679 = 0.1e1 / t703;
t697 = sin(qJ(2,1));
t818 = t679 * t697;
t687 = pkin(3) + qJ(3,3);
t656 = 0.1e1 / t687;
t688 = pkin(3) + qJ(3,2);
t659 = 0.1e1 / t688;
t689 = pkin(3) + qJ(3,1);
t662 = 0.1e1 / t689;
t708 = pkin(1) + pkin(2);
t684 = 1 / t708;
t657 = 0.1e1 / t687 ^ 2;
t660 = 0.1e1 / t688 ^ 2;
t663 = 0.1e1 / t689 ^ 2;
t848 = pkin(1) * t684;
t694 = sin(qJ(1,3));
t700 = cos(qJ(1,3));
t805 = t699 * t708;
t644 = t694 * t687 + t700 * t805;
t690 = legFrame(3,2);
t650 = sin(t690);
t653 = cos(t690);
t705 = xDP(3);
t706 = xDP(2);
t707 = xDP(1);
t812 = t694 * t699;
t817 = t687 * t700;
t617 = (-t707 * t817 + (t693 * t706 + t707 * t812) * t708) * t653 + (t706 * t817 - (-t693 * t707 + t706 * t812) * t708) * t650 + t644 * t705;
t641 = t650 * t707 + t653 * t706;
t814 = t693 * t641;
t623 = (t700 * t705 + (-t650 * t706 + t653 * t707) * t694) * t699 + t814;
t638 = t641 ^ 2;
t658 = t656 * t657;
t668 = t699 ^ 2;
t671 = t669 / t668;
t827 = t657 * t669;
t593 = -t623 * t658 * t669 * t617 + (-(-t623 * t708 + t617) * t623 * t827 + t684 * t638 * t671) * t656;
t847 = t593 * qJ(3,3);
t696 = sin(qJ(1,2));
t702 = cos(qJ(1,2));
t804 = t701 * t708;
t645 = t696 * t688 + t702 * t804;
t691 = legFrame(2,2);
t651 = sin(t691);
t654 = cos(t691);
t809 = t696 * t701;
t816 = t688 * t702;
t618 = (-t707 * t816 + (t695 * t706 + t707 * t809) * t708) * t654 + (t706 * t816 - (-t695 * t707 + t706 * t809) * t708) * t651 + t645 * t705;
t642 = t651 * t707 + t654 * t706;
t811 = t695 * t642;
t624 = (t702 * t705 + (-t651 * t706 + t654 * t707) * t696) * t701 + t811;
t639 = t642 ^ 2;
t661 = t659 * t660;
t673 = t701 ^ 2;
t676 = t674 / t673;
t824 = t660 * t674;
t594 = -t624 * t661 * t674 * t618 + (-(-t624 * t708 + t618) * t624 * t824 + t684 * t639 * t676) * t659;
t846 = t594 * qJ(3,2);
t698 = sin(qJ(1,1));
t704 = cos(qJ(1,1));
t803 = t703 * t708;
t646 = t698 * t689 + t704 * t803;
t692 = legFrame(1,2);
t652 = sin(t692);
t655 = cos(t692);
t806 = t698 * t703;
t815 = t689 * t704;
t619 = (-t707 * t815 + (t697 * t706 + t707 * t806) * t708) * t655 + (t706 * t815 - (-t697 * t707 + t706 * t806) * t708) * t652 + t646 * t705;
t643 = t652 * t707 + t655 * t706;
t808 = t697 * t643;
t625 = (t704 * t705 + (-t652 * t706 + t655 * t707) * t698) * t703 + t808;
t640 = t643 ^ 2;
t664 = t662 * t663;
t678 = t703 ^ 2;
t681 = t679 / t678;
t821 = t663 * t679;
t595 = -t625 * t664 * t679 * t619 + (-(-t625 * t708 + t619) * t625 * t821 + t684 * t640 * t681) * t662;
t845 = t595 * qJ(3,1);
t670 = 0.1e1 / t699 ^ 2;
t620 = t623 ^ 2;
t838 = t620 * t657;
t844 = (t670 - 0.2e1) * t838 * t669;
t675 = 0.1e1 / t701 ^ 2;
t621 = t624 ^ 2;
t837 = t621 * t660;
t843 = (t675 - 0.2e1) * t837 * t674;
t680 = 0.1e1 / t703 ^ 2;
t622 = t625 ^ 2;
t836 = t622 * t663;
t842 = (t680 - 0.2e1) * t836 * t679;
t841 = t617 * t708;
t840 = t618 * t708;
t839 = t619 * t708;
t835 = t623 * t657;
t834 = t624 * t660;
t833 = t625 * t663;
t685 = 1 / t708 ^ 2;
t832 = t638 * t685;
t831 = t639 * t685;
t830 = t640 * t685;
t829 = t656 * t693;
t828 = t656 * t700;
t826 = t659 * t695;
t825 = t659 * t702;
t823 = t662 * t697;
t822 = t662 * t704;
t813 = t693 * t708;
t810 = t695 * t708;
t807 = t697 * t708;
t683 = t708 ^ 2;
t727 = (-t623 * t693 + t641) * t669 ^ 2 * t656 * t641 + (-((-t623 + t814) * t687 * t669 + (-t623 * t683 + t841) * t699 * t656) * t827 - t658 * t841) * t623;
t733 = -pkin(1) * t832 + 0.2e1 * t617 * t835;
t783 = t623 * t829;
t745 = t641 * t684 * t783;
t665 = t693 ^ 2;
t778 = t671 * t832;
t754 = t665 * t778;
t802 = (qJ(3,3) ^ 2 + t668 * t709) * t593 + (-qJ(3,3) * t754 - t699 * t727) * pkin(1) + (t733 * qJ(3,3) + t745 * t849) * t669;
t726 = (-t624 * t695 + t642) * t674 ^ 2 * t659 * t642 + (-((-t624 + t811) * t688 * t674 + (-t624 * t683 + t840) * t701 * t659) * t824 - t661 * t840) * t624;
t732 = -pkin(1) * t831 + 0.2e1 * t618 * t834;
t781 = t624 * t826;
t744 = t642 * t684 * t781;
t666 = t695 ^ 2;
t776 = t676 * t831;
t753 = t666 * t776;
t801 = (qJ(3,2) ^ 2 + t673 * t709) * t594 + (-qJ(3,2) * t753 - t701 * t726) * pkin(1) + (t732 * qJ(3,2) + t744 * t849) * t674;
t725 = (-t625 * t697 + t643) * t679 ^ 2 * t662 * t643 + (-((-t625 + t808) * t689 * t679 + (-t625 * t683 + t839) * t703 * t662) * t821 - t664 * t839) * t625;
t731 = -pkin(1) * t830 + 0.2e1 * t619 * t833;
t779 = t625 * t823;
t743 = t643 * t684 * t779;
t667 = t697 ^ 2;
t774 = t681 * t830;
t752 = t667 * t774;
t800 = (qJ(3,1) ^ 2 + t678 * t709) * t595 + (-qJ(3,1) * t752 - t703 * t725) * pkin(1) + (t731 * qJ(3,1) + t743 * t849) * t679;
t799 = -t699 * t593 * pkin(1) + (-qJ(3,3) * t838 + t745 * t850) * t670 + t727;
t798 = -t701 * t594 * pkin(1) + (-qJ(3,2) * t837 + t744 * t850) * t675 + t726;
t797 = -t703 * t595 * pkin(1) + (-qJ(3,1) * t836 + t743 * t850) * t680 + t725;
t796 = 0.2e1 * t623 * t641;
t795 = 0.2e1 * t624 * t642;
t794 = 0.2e1 * t625 * t643;
t793 = t593 * t656 * t669;
t792 = t593 * t820;
t791 = t594 * t659 * t674;
t790 = t594 * t819;
t789 = t595 * t662 * t679;
t788 = t595 * t818;
t787 = t620 * t658 * t670;
t786 = t621 * t661 * t675;
t785 = t622 * t664 * t680;
t784 = t641 * t835;
t782 = t642 * t834;
t780 = t643 * t833;
t672 = 0.1e1 / t668 ^ 2;
t777 = t638 * t672 * t693;
t677 = 0.1e1 / t673 ^ 2;
t775 = t639 * t677 * t695;
t682 = 0.1e1 / t678 ^ 2;
t773 = t640 * t682 * t697;
t772 = t657 * t670 * t693;
t771 = t660 * t675 * t695;
t770 = t663 * t680 * t697;
t769 = t593 * t828;
t768 = t594 * t825;
t767 = t595 * t822;
t766 = t656 * (-pkin(1) * t754 + t733 * t669 + 0.2e1 * t847);
t765 = t659 * (-pkin(1) * t753 + t732 * t674 + 0.2e1 * t846);
t764 = t662 * (-pkin(1) * t752 + t731 * t679 + 0.2e1 * t845);
t763 = t802 * t669;
t762 = t801 * t674;
t761 = t800 * t679;
t760 = 0.2e1 * t593 * t829;
t759 = 0.2e1 * t594 * t826;
t758 = 0.2e1 * t595 * t823;
t757 = t700 * t784;
t756 = t702 * t782;
t755 = t704 * t780;
t751 = t665 * t793;
t750 = t666 * t791;
t749 = t667 * t789;
t748 = t669 * t766;
t747 = t674 * t765;
t746 = t679 * t764;
t647 = 0.2e1 * t668 - 0.1e1;
t742 = 0.2e1 * t647 * t671 * t784;
t648 = 0.2e1 * t673 - 0.1e1;
t741 = 0.2e1 * t648 * t676 * t782;
t649 = 0.2e1 * t678 - 0.1e1;
t740 = 0.2e1 * t649 * t681 * t780;
t739 = (t665 * t672 + t670) * t656 * t638;
t738 = (t666 * t677 + t675) * t659 * t639;
t737 = (t667 * t682 + t680) * t662 * t640;
t736 = t694 * t805 - t817;
t735 = t696 * t804 - t816;
t734 = t698 * t803 - t815;
t730 = (pkin(1) * t623 - 0.2e1 * t617) * t656 * t670 * t783 + (pkin(1) * t778 - t847) * t820;
t729 = (pkin(1) * t624 - 0.2e1 * t618) * t659 * t675 * t781 + (pkin(1) * t776 - t846) * t819;
t728 = (pkin(1) * t625 - 0.2e1 * t619) * t662 * t680 * t779 + (pkin(1) * t774 - t845) * t818;
t724 = t650 * t792 + t651 * t790 + t652 * t788;
t723 = t653 * t792 + t654 * t790 + t655 * t788;
t686 = t684 / t683;
t637 = t697 * t652 + t655 * t806;
t636 = t695 * t651 + t654 * t809;
t635 = t693 * t650 + t653 * t812;
t634 = -t652 * t806 + t697 * t655;
t633 = -t651 * t809 + t695 * t654;
t632 = -t650 * t812 + t693 * t653;
t631 = t652 * t807 + t734 * t655;
t630 = t651 * t810 + t735 * t654;
t629 = t650 * t813 + t736 * t653;
t628 = -t735 * t651 + t654 * t810;
t627 = -t734 * t652 + t655 * t807;
t626 = -t736 * t650 + t653 * t813;
t1 = [t635 * t793 + t636 * t791 + t637 * t789, 0, 0, t635 * t751 + t636 * t750 + t637 * t749 + ((-t622 * t652 + t637 * t794) * t770 + (-t621 * t651 + t636 * t795) * t771 + (-t620 * t650 + t635 * t796) * t772) * t684, t635 * t760 + t636 * t759 + t637 * t758 + (t635 * t742 + t636 * t741 + t637 * t740 + t650 * t844 + t651 * t843 + t652 * t842) * t684, t724 * t684 + (t635 * t739 + t636 * t738 + t637 * t737) * t685, (t593 * t650 + t594 * t651 + t595 * t652) * t684, (t650 * t777 + t651 * t775 + t652 * t773) * t686, 0, 0, -t629 * t787 - t630 * t786 - t631 * t785 + t635 * t748 + t636 * t747 + t637 * t746 - t724 * t848, (t631 * t797 + t637 * t761) * t662 + (t630 * t798 + t636 * t762) * t659 + (t629 * t799 + t635 * t763) * t656 + (t650 * t730 + t651 * t729 + t652 * t728) * t848, 0; t632 * t793 + t633 * t791 + t634 * t789, 0, 0, t632 * t751 + t633 * t750 + t634 * t749 + ((-t622 * t655 + t634 * t794) * t770 + (-t621 * t654 + t633 * t795) * t771 + (-t620 * t653 + t632 * t796) * t772) * t684, t632 * t760 + t633 * t759 + t634 * t758 + (t632 * t742 + t633 * t741 + t634 * t740 + t653 * t844 + t654 * t843 + t655 * t842) * t684, t723 * t684 + (t632 * t739 + t633 * t738 + t634 * t737) * t685, (t593 * t653 + t594 * t654 + t595 * t655) * t684, (t653 * t777 + t654 * t775 + t655 * t773) * t686, 0, 0, -t626 * t787 - t627 * t785 - t628 * t786 + t632 * t748 + t633 * t747 + t634 * t746 - t723 * t848, (t627 * t797 + t634 * t761) * t662 + (t628 * t798 + t633 * t762) * t659 + (t626 * t799 + t632 * t763) * t656 + (t653 * t730 + t654 * t729 + t655 * t728) * t848, 0; t767 + t768 + t769, 0, 0, t665 * t769 + t666 * t768 + t667 * t767 + 0.2e1 * (t755 * t818 + t756 * t819 + t757 * t820) * t684, 0.2e1 * t693 * t699 * t769 + 0.2e1 * t695 * t701 * t768 + 0.2e1 * t697 * t703 * t767 + 0.2e1 * (t647 * t670 * t757 + t648 * t675 * t756 + t649 * t680 * t755) * t684, ((t667 * t681 + t679) * t640 * t822 + (t666 * t676 + t674) * t639 * t825 + (t665 * t671 + t669) * t638 * t828) * t685, 0, 0, 0, 0, -t644 * t787 - t645 * t786 - t646 * t785 + t700 * t766 + t702 * t765 + t704 * t764, (t646 * t797 + t704 * t800) * t662 + (t645 * t798 + t702 * t801) * t659 + (t644 * t799 + t700 * t802) * t656, 0;];
tau_reg  = t1;
