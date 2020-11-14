% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:35
% EndTime: 2020-03-09 21:16:37
% DurationCPUTime: 2.08s
% Computational Cost: add. (1718->209), mult. (5637->517), div. (2193->17), fcn. (5130->18), ass. (0->249)
t673 = cos(qJ(3,3));
t646 = 0.1e1 / t673;
t695 = t673 ^ 2;
t647 = 0.1e1 / t695;
t648 = t646 * t647;
t681 = 0.1e1 / pkin(2);
t664 = legFrame(3,2);
t625 = sin(t664);
t628 = cos(t664);
t680 = xDP(1);
t864 = xDP(2);
t613 = -t625 * t864 + t628 * t680;
t610 = t613 ^ 2;
t668 = sin(qJ(2,3));
t633 = 0.1e1 / t668;
t851 = t610 * t633;
t679 = xDP(3);
t667 = sin(qJ(3,3));
t674 = cos(qJ(2,3));
t819 = t667 * t674;
t601 = t613 * t819 - t673 * t679;
t598 = t601 ^ 2;
t634 = 0.1e1 / t668 ^ 2;
t860 = t598 * t634;
t859 = t633 * t860;
t589 = (t851 + t859) * t681 * t648;
t632 = t667 ^ 2;
t651 = t674 ^ 2;
t841 = t633 * t651;
t766 = t647 * t841;
t880 = t589 * (t632 * t766 - t668);
t879 = t589 * (t668 + t841);
t675 = cos(qJ(3,2));
t652 = 0.1e1 / t675;
t699 = t675 ^ 2;
t653 = 0.1e1 / t699;
t654 = t652 * t653;
t665 = legFrame(2,2);
t626 = sin(t665);
t629 = cos(t665);
t615 = -t626 * t864 + t629 * t680;
t611 = t615 ^ 2;
t670 = sin(qJ(2,2));
t638 = 0.1e1 / t670;
t849 = t611 * t638;
t669 = sin(qJ(3,2));
t676 = cos(qJ(2,2));
t815 = t669 * t676;
t602 = t615 * t815 - t675 * t679;
t599 = t602 ^ 2;
t639 = 0.1e1 / t670 ^ 2;
t858 = t599 * t639;
t857 = t638 * t858;
t590 = (t849 + t857) * t681 * t654;
t637 = t669 ^ 2;
t657 = t676 ^ 2;
t837 = t638 * t657;
t765 = t653 * t837;
t878 = t590 * (t637 * t765 - t670);
t877 = t590 * (t670 + t837);
t677 = cos(qJ(3,1));
t658 = 0.1e1 / t677;
t703 = t677 ^ 2;
t659 = 0.1e1 / t703;
t660 = t658 * t659;
t666 = legFrame(1,2);
t627 = sin(t666);
t630 = cos(t666);
t617 = -t627 * t864 + t630 * t680;
t612 = t617 ^ 2;
t672 = sin(qJ(2,1));
t643 = 0.1e1 / t672;
t847 = t612 * t643;
t671 = sin(qJ(3,1));
t678 = cos(qJ(2,1));
t811 = t671 * t678;
t603 = t617 * t811 - t677 * t679;
t600 = t603 ^ 2;
t644 = 0.1e1 / t672 ^ 2;
t856 = t600 * t644;
t855 = t643 * t856;
t591 = (t847 + t855) * t681 * t660;
t642 = t671 ^ 2;
t663 = t678 ^ 2;
t833 = t643 * t663;
t764 = t659 * t833;
t876 = t591 * (t642 * t764 - t672);
t875 = t591 * (t672 + t833);
t649 = 0.1e1 / t695 ^ 2;
t682 = 0.1e1 / pkin(2) ^ 2;
t789 = t682 * t860;
t595 = t649 * t789;
t820 = t667 * t668;
t840 = t633 * t674;
t845 = t613 * t674;
t854 = t601 * t634;
t865 = t681 ^ 2;
t577 = (-(t601 * t840 + t613 * t820) * t854 + (-t601 * t667 - t845) * t633 * t613) * t648 * t646 * t865;
t809 = t674 * t577;
t798 = t633 * t809;
t874 = t646 * (0.2e1 * t632 * t798 + 0.2e1 * t647 * t789 - t595);
t655 = 0.1e1 / t699 ^ 2;
t787 = t682 * t858;
t596 = t655 * t787;
t816 = t669 * t670;
t836 = t638 * t676;
t844 = t615 * t676;
t853 = t602 * t639;
t578 = (-(t602 * t836 + t615 * t816) * t853 + (-t602 * t669 - t844) * t638 * t615) * t654 * t652 * t865;
t808 = t676 * t578;
t796 = t638 * t808;
t873 = t652 * (0.2e1 * t637 * t796 + 0.2e1 * t653 * t787 - t596);
t661 = 0.1e1 / t703 ^ 2;
t785 = t682 * t856;
t597 = t661 * t785;
t812 = t671 * t672;
t832 = t643 * t678;
t843 = t617 * t678;
t852 = t603 * t644;
t579 = (-(t603 * t832 + t617 * t812) * t852 + (-t603 * t671 - t843) * t643 * t617) * t660 * t658 * t865;
t807 = t678 * t579;
t794 = t643 * t807;
t872 = t658 * (0.2e1 * t642 * t794 + 0.2e1 * t659 * t785 - t597);
t758 = 0.2e1 * t603 * t843;
t831 = t644 * t661;
t871 = (-t600 * t671 + t642 * t758) * t831;
t759 = 0.2e1 * t602 * t844;
t835 = t639 * t655;
t870 = (-t599 * t669 + t637 * t759) * t835;
t760 = 0.2e1 * t601 * t845;
t839 = t634 * t649;
t869 = (-t598 * t667 + t632 * t760) * t839;
t641 = t671 * t642;
t662 = t658 * t661;
t813 = t671 * t660;
t868 = (t641 * t662 + t813) * t612 * t832;
t636 = t669 * t637;
t656 = t652 * t655;
t817 = t669 * t654;
t867 = (t636 * t656 + t817) * t611 * t836;
t631 = t667 * t632;
t650 = t646 * t649;
t821 = t667 * t648;
t866 = (t631 * t650 + t821) * t610 * t840;
t863 = t589 * t681;
t862 = t590 * t681;
t861 = t591 * t681;
t850 = t610 * t682;
t848 = t611 * t682;
t846 = t612 * t682;
t842 = t633 * t646;
t838 = t638 * t652;
t834 = t643 * t658;
t830 = t646 * t577;
t829 = t646 * t667;
t828 = t650 * t674;
t827 = t652 * t578;
t826 = t652 * t669;
t825 = t656 * t676;
t824 = t658 * t579;
t823 = t658 * t671;
t822 = t662 * t678;
t818 = t668 * t673;
t814 = t670 * t675;
t810 = t672 * t677;
t779 = t647 * t850;
t592 = t595 + t779;
t721 = t633 * t682 * t760;
t806 = -t632 * t668 * t648 * t850 - t592 * t818 + t673 * t809 + t721 * t821;
t776 = t653 * t848;
t593 = t596 + t776;
t720 = t638 * t682 * t759;
t805 = -t637 * t670 * t654 * t848 - t593 * t814 + t675 * t808 + t720 * t817;
t773 = t659 * t846;
t594 = t597 + t773;
t719 = t643 * t682 * t758;
t804 = -t642 * t672 * t660 * t846 - t594 * t810 + t677 * t807 + t719 * t813;
t803 = (-t668 * t779 - t809) * t667 + t592 * t820 + t647 * t721;
t802 = (-t670 * t776 - t808) * t669 + t593 * t816 + t653 * t720;
t801 = (-t672 * t773 - t807) * t671 + t594 * t812 + t659 * t719;
t683 = t681 * t682;
t800 = 0.2e1 * t683;
t799 = t633 * t830;
t797 = t638 * t827;
t795 = t643 * t824;
t793 = t589 * t842;
t792 = t590 * t838;
t791 = t591 * t834;
t790 = t650 * t860;
t788 = t656 * t858;
t786 = t662 * t856;
t784 = t613 * t854;
t783 = t615 * t853;
t782 = t617 * t852;
t780 = t610 * t649 * t667;
t777 = t611 * t655 * t669;
t774 = t612 * t661 * t671;
t772 = t625 * t829;
t771 = t626 * t826;
t770 = t627 * t823;
t769 = t628 * t829;
t768 = t629 * t826;
t767 = t630 * t823;
t763 = t589 * t840;
t762 = t590 * t836;
t761 = t591 * t832;
t754 = t647 * t798;
t753 = t653 * t796;
t752 = t659 * t794;
t751 = t589 * t647 * t819;
t750 = t590 * t653 * t815;
t749 = t591 * t659 * t811;
t748 = t828 * t859;
t747 = t825 * t857;
t746 = t822 * t855;
t745 = t667 * t784;
t744 = t669 * t783;
t743 = t671 * t782;
t742 = t646 * t798;
t741 = t652 * t796;
t740 = t658 * t794;
t739 = t806 * t842;
t738 = t805 * t838;
t737 = t804 * t834;
t736 = t803 * t842;
t735 = t802 * t838;
t734 = t801 * t834;
t733 = t631 * t754;
t732 = t667 * t754;
t731 = t636 * t753;
t730 = t669 * t753;
t729 = t641 * t752;
t728 = t671 * t752;
t727 = t589 * t667 * t766;
t726 = t590 * t669 * t765;
t725 = t591 * t671 * t764;
t622 = -0.1e1 + 0.2e1 * t695;
t712 = t622 * t745 * t828;
t623 = -0.1e1 + 0.2e1 * t699;
t711 = t623 * t744 * t825;
t624 = -0.1e1 + 0.2e1 * t703;
t710 = t624 * t743 * t822;
t609 = t627 * t671 + t630 * t810;
t608 = t626 * t669 + t629 * t814;
t607 = t625 * t667 + t628 * t818;
t606 = t627 * t810 - t630 * t671;
t605 = t626 * t814 - t629 * t669;
t604 = t625 * t818 - t628 * t667;
t1 = [t604 * t793 + t605 * t792 + t606 * t791, (t628 * t732 + t629 * t730 + t630 * t728) * t681, t604 * t742 + t605 * t741 + t606 * t740 + (-t604 * t790 - t605 * t788 - t606 * t786) * t682 + (t628 * t727 + t629 * t726 + t630 * t725) * t681, -t604 * t830 - t605 * t827 - t606 * t824 + (-t604 * t748 - t605 * t747 - t606 * t746) * t682 + (-t628 * t751 - t629 * t750 - t630 * t749) * t681, (t628 * t733 + t629 * t731 + t630 * t729) * t681 + (-t628 * t869 - t629 * t870 - t630 * t871) * t683, (-t628 * t712 - t629 * t711 - t630 * t710) * t800 + (t628 * t874 + t629 * t873 + t630 * t872) * t681, (-t577 * t769 - t578 * t768 - t579 * t767) * t681 + (t628 * t866 + t629 * t867 + t630 * t868) * t683, (-t577 * t628 - t578 * t629 - t579 * t630) * t681, (-t628 * t780 - t629 * t777 - t630 * t774) * t683, t606 * t737 + t605 * t738 + t604 * t739 + (t767 * t875 + t768 * t877 + t769 * t879) * t681, t606 * t734 + t605 * t735 + t604 * t736 + (-t628 * t880 - t629 * t878 - t630 * t876) * t681, 0; t607 * t793 + t608 * t792 + t609 * t791, (-t625 * t732 - t626 * t730 - t627 * t728) * t681, t607 * t742 + t608 * t741 + t609 * t740 + (-t607 * t790 - t608 * t788 - t609 * t786) * t682 + (-t625 * t727 - t626 * t726 - t627 * t725) * t681, -t607 * t830 - t608 * t827 - t609 * t824 + (-t607 * t748 - t608 * t747 - t609 * t746) * t682 + (t625 * t751 + t626 * t750 + t627 * t749) * t681, (-t625 * t733 - t626 * t731 - t627 * t729) * t681 + (t625 * t869 + t626 * t870 + t627 * t871) * t683, (t625 * t712 + t626 * t711 + t627 * t710) * t800 + (-t625 * t874 - t626 * t873 - t627 * t872) * t681, (t577 * t772 + t578 * t771 + t579 * t770) * t681 + (-t625 * t866 - t626 * t867 - t627 * t868) * t683, (t577 * t625 + t578 * t626 + t579 * t627) * t681, (t625 * t780 + t626 * t777 + t627 * t774) * t683, t609 * t737 + t608 * t738 + t607 * t739 + (-t770 * t875 - t771 * t877 - t772 * t879) * t681, t609 * t734 + t608 * t735 + t607 * t736 + (t625 * t880 + t626 * t878 + t627 * t876) * t681, 0; t761 + t762 + t763, (-t795 - t797 - t799) * t681, t577 * t841 + t578 * t837 + t579 * t833 + (-t598 * t674 * t839 - t599 * t676 * t835 - t600 * t678 * t831) * t682 + (-t646 * t763 - t652 * t762 - t658 * t761) * t681, -t809 - t808 - t807 + (-t649 * t651 * t859 - t655 * t657 * t857 - t661 * t663 * t855) * t682 + (t589 * t646 + t590 * t652 + t591 * t658) * t681, (t648 * t745 + t654 * t744 + t660 * t743) * t800 + (-t632 * t799 - t637 * t797 - t642 * t795) * t681, 0.2e1 * (t622 * t649 * t784 + t623 * t655 * t783 + t624 * t661 * t782) * t683 + 0.2e1 * (-t577 * t633 * t667 - t578 * t638 * t669 - t579 * t643 * t671) * t681, ((-t642 * t661 - t659) * t847 + (-t637 * t655 - t653) * t849 + (-t632 * t649 - t647) * t851) * t683, 0, 0, (t804 - t861) * t832 + (t805 - t862) * t836 + (t806 - t863) * t840, (t823 * t861 + t801) * t832 + (t826 * t862 + t802) * t836 + (t829 * t863 + t803) * t840, 0;];
tau_reg  = t1;
