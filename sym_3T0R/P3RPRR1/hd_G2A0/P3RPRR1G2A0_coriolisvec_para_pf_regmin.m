% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRR1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:25:05
% EndTime: 2020-03-09 21:25:07
% DurationCPUTime: 2.74s
% Computational Cost: add. (34308->233), mult. (26523->357), div. (1473->10), fcn. (13608->59), ass. (0->203)
t778 = pkin(7) + qJ(3,3);
t786 = sin(qJ(3,3));
t704 = pkin(1) * sin(t778) + t786 * pkin(2);
t696 = 0.1e1 / t704 ^ 2;
t783 = legFrame(3,2);
t847 = qJ(1,3) + pkin(7);
t740 = t783 + t847;
t734 = qJ(3,3) + t740;
t741 = -t783 + t847;
t735 = qJ(3,3) + t741;
t674 = sin(t734) + sin(t735);
t677 = -cos(t735) + cos(t734);
t752 = cos(qJ(1,3) + t778);
t799 = xDP(2);
t800 = xDP(1);
t920 = 2 * xDP(3);
t659 = t674 * t800 + t677 * t799 + t752 * t920;
t872 = t659 / 0.2e1;
t895 = t696 * t872;
t779 = pkin(7) + qJ(3,2);
t787 = sin(qJ(3,2));
t705 = pkin(1) * sin(t779) + t787 * pkin(2);
t699 = 0.1e1 / t705 ^ 2;
t784 = legFrame(2,2);
t848 = qJ(1,2) + pkin(7);
t742 = t784 + t848;
t736 = qJ(3,2) + t742;
t743 = -t784 + t848;
t737 = qJ(3,2) + t743;
t675 = sin(t736) + sin(t737);
t678 = -cos(t737) + cos(t736);
t753 = cos(qJ(1,2) + t779);
t661 = t675 * t800 + t678 * t799 + t753 * t920;
t868 = t661 / 0.2e1;
t897 = t699 * t868;
t780 = pkin(7) + qJ(3,1);
t788 = sin(qJ(3,1));
t706 = pkin(1) * sin(t780) + t788 * pkin(2);
t702 = 0.1e1 / t706 ^ 2;
t785 = legFrame(1,2);
t849 = qJ(1,1) + pkin(7);
t744 = t785 + t849;
t738 = qJ(3,1) + t744;
t745 = -t785 + t849;
t739 = qJ(3,1) + t745;
t676 = sin(t738) + sin(t739);
t679 = -cos(t739) + cos(t738);
t754 = cos(qJ(1,1) + t780);
t660 = t676 * t800 + t679 * t799 + t754 * t920;
t870 = t660 / 0.2e1;
t896 = t702 * t870;
t771 = qJ(1,3) + t783;
t772 = qJ(1,3) - t783;
t662 = -t674 * pkin(3) + (-sin(t740) - sin(t741)) * pkin(2) + (-sin(t771) - sin(t772)) * pkin(1);
t665 = -t677 * pkin(3) + (-cos(t740) + cos(t741)) * pkin(2) + (-cos(t771) + cos(t772)) * pkin(1);
t671 = -cos(qJ(1,3)) * pkin(1) - pkin(2) * cos(t847) - pkin(3) * t752;
t695 = 0.1e1 / t704;
t802 = 0.1e1 / pkin(3);
t878 = (t662 * t800 + t665 * t799 + t671 * t920) * t695 * t802;
t648 = t878 / 0.2e1;
t653 = t695 * t872;
t638 = t653 + t648 / 0.2e1;
t782 = cos(pkin(7));
t761 = t782 * pkin(1) + pkin(2);
t789 = cos(qJ(3,3));
t781 = sin(pkin(7));
t892 = pkin(1) * t781;
t689 = t786 * t761 + t789 * t892;
t837 = t695 * t878;
t823 = t638 * t689 * t837;
t644 = t653 + t648;
t890 = pkin(3) * t644;
t633 = t695 * t648 * t890;
t765 = cos(t778);
t891 = pkin(2) * t789;
t635 = -t890 + (-pkin(1) * t765 - t891) * t653;
t855 = t781 * t786;
t668 = -t891 + (-t782 * t789 + t855) * pkin(1);
t803 = pkin(1) ^ 2;
t777 = pkin(2) ^ 2 + t803;
t801 = pkin(3) ^ 2;
t850 = 0.2e1 * pkin(1);
t851 = 0.2e1 * pkin(2) * pkin(3);
t894 = pkin(2) * t782;
t882 = (t638 * t789 * t851 + t777 * t653 + t644 * t801 + (pkin(3) * t638 * t765 + t653 * t894) * t850) * t802;
t686 = -pkin(1) * t855 + t789 * t761;
t817 = t644 * (pkin(3) + t686) / t689 * t878;
t913 = -t817 / 0.2e1;
t829 = t695 * (0.2e1 * t633 + t913 + (-0.2e1 * t635 - t882) * t895) * t668;
t932 = -t823 / 0.2e1 - t829 / 0.2e1;
t775 = qJ(1,1) + t785;
t776 = qJ(1,1) - t785;
t664 = -t676 * pkin(3) + (-sin(t744) - sin(t745)) * pkin(2) + (-sin(t775) - sin(t776)) * pkin(1);
t667 = -t679 * pkin(3) + (-cos(t744) + cos(t745)) * pkin(2) + (-cos(t775) + cos(t776)) * pkin(1);
t673 = -cos(qJ(1,1)) * pkin(1) - pkin(2) * cos(t849) - pkin(3) * t754;
t701 = 0.1e1 / t706;
t877 = (t664 * t800 + t667 * t799 + t673 * t920) * t701 * t802;
t649 = t877 / 0.2e1;
t654 = t701 * t870;
t639 = t654 + t649 / 0.2e1;
t793 = cos(qJ(3,1));
t691 = t788 * t761 + t793 * t892;
t836 = t701 * t877;
t822 = t639 * t691 * t836;
t645 = t654 + t649;
t888 = t645 * pkin(3);
t632 = t701 * t649 * t888;
t767 = cos(t780);
t886 = t793 * pkin(2);
t637 = -t888 + (-pkin(1) * t767 - t886) * t654;
t853 = t781 * t788;
t670 = t886 + (t782 * t793 - t853) * pkin(1);
t881 = (t639 * t793 * t851 + t777 * t654 + t645 * t801 + (pkin(3) * t639 * t767 + t654 * t894) * t850) * t802;
t688 = -pkin(1) * t853 + t793 * t761;
t816 = t645 * (pkin(3) + t688) / t691 * t877;
t914 = -t816 / 0.2e1;
t827 = t701 * (0.2e1 * t632 + t914 + (-0.2e1 * t637 - t881) * t896) * t670;
t931 = -t822 / 0.2e1 + t827 / 0.2e1;
t773 = qJ(1,2) + t784;
t774 = qJ(1,2) - t784;
t663 = -t675 * pkin(3) + (-sin(t742) - sin(t743)) * pkin(2) + (-sin(t773) - sin(t774)) * pkin(1);
t666 = -t678 * pkin(3) + (-cos(t742) + cos(t743)) * pkin(2) + (-cos(t773) + cos(t774)) * pkin(1);
t672 = -cos(qJ(1,2)) * pkin(1) - pkin(2) * cos(t848) - pkin(3) * t753;
t698 = 0.1e1 / t705;
t879 = (t663 * t800 + t666 * t799 + t672 * t920) * t698 * t802;
t647 = t879 / 0.2e1;
t655 = t698 * t868;
t640 = t655 + t647 / 0.2e1;
t791 = cos(qJ(3,2));
t690 = t787 * t761 + t791 * t892;
t838 = t698 * t879;
t821 = t640 * t690 * t838;
t646 = t655 + t647;
t889 = pkin(3) * t646;
t634 = t698 * t647 * t889;
t766 = cos(t779);
t887 = t791 * pkin(2);
t636 = -t889 + (-pkin(1) * t766 - t887) * t655;
t854 = t781 * t787;
t669 = -t887 + (-t782 * t791 + t854) * pkin(1);
t880 = (t640 * t791 * t851 + t777 * t655 + t646 * t801 + (pkin(3) * t640 * t766 + t655 * t894) * t850) * t802;
t687 = -pkin(1) * t854 + t791 * t761;
t815 = t646 * (pkin(3) + t687) / t690 * t879;
t915 = -t815 / 0.2e1;
t828 = t698 * (0.2e1 * t634 + t915 + (-0.2e1 * t636 - t880) * t897) * t669;
t930 = -t821 / 0.2e1 - t828 / 0.2e1;
t929 = t695 / 0.2e1;
t928 = t698 / 0.2e1;
t927 = t701 / 0.2e1;
t820 = (t647 + 0.2e1 * t655) * t669 * t838;
t860 = t690 * t698;
t825 = -0.2e1 * (t634 - t815 / 0.4e1 + (-t636 - t880 / 0.2e1) * t897) * t860;
t923 = t820 / 0.4e1 + t825 / 0.2e1;
t819 = (t648 + 0.2e1 * t653) * t668 * t837;
t861 = t689 * t695;
t826 = -0.2e1 * (t633 - t817 / 0.4e1 + (-t635 - t882 / 0.2e1) * t895) * t861;
t922 = t819 / 0.4e1 + t826 / 0.2e1;
t818 = (t649 + 0.2e1 * t654) * t670 * t836;
t859 = t691 * t701;
t824 = -0.2e1 * (t632 - t816 / 0.4e1 + (-t637 - t881 / 0.2e1) * t896) * t859;
t921 = -t818 / 0.4e1 + t824 / 0.2e1;
t919 = t671 / 0.4e1;
t918 = t672 / 0.4e1;
t917 = t673 / 0.4e1;
t916 = t802 / 0.2e1;
t903 = t679 * t927;
t902 = t678 * t928;
t901 = t677 * t929;
t900 = t676 * t927;
t899 = t675 * t928;
t898 = t674 * t929;
t620 = t633 + t913 + (-t635 - t882) * t895;
t885 = t620 * t695;
t622 = t634 + t915 + (-t636 - t880) * t897;
t884 = t622 * t698;
t624 = t632 + t914 + (-t637 - t881) * t896;
t883 = t624 * t701;
t876 = t659 ^ 2 * t695 * t696;
t875 = t660 ^ 2 * t701 * t702;
t874 = t661 ^ 2 * t698 * t699;
t858 = t695 * t752;
t857 = t698 * t753;
t856 = t701 * t754;
t852 = t802 / 0.8e1;
t629 = -t635 * t895 + t633;
t844 = t629 * t686 * t695;
t843 = t629 * t861;
t630 = -t636 * t897 + t634;
t842 = t630 * t687 * t698;
t841 = t630 * t860;
t631 = -t637 * t896 + t632;
t840 = t631 * t688 * t701;
t839 = t631 * t859;
t835 = t686 * t876;
t834 = t689 * t876;
t833 = t688 * t875;
t832 = t691 * t875;
t831 = t687 * t874;
t830 = t690 * t874;
t814 = t629 * t858 + t630 * t857 + t631 * t856;
t813 = t629 * t898 + t630 * t899 + t631 * t900;
t812 = t629 * t901 + t630 * t902 + t631 * t903;
t1 = [t813, 0, 0, t803 * t813, t620 * t898 + t622 * t899 + t624 * t900 + (t662 * t885 + t663 * t884 + t664 * t883) * t916, (t662 * t834 + t663 * t830 + t664 * t832) * t852 + (t662 * t844 + t663 * t842 + t664 * t840) * t916 + t931 * t676 + t930 * t675 + t932 * t674, (t662 * t835 + t663 * t831 + t664 * t833) * t852 + (-t662 * t843 - t663 * t841 - t664 * t839) * t916 + t921 * t676 + t923 * t675 + t922 * t674, 0; t812, 0, 0, t803 * t812, t620 * t901 + t622 * t902 + t624 * t903 + (t665 * t885 + t666 * t884 + t667 * t883) * t916, (t665 * t834 + t666 * t830 + t667 * t832) * t852 + (t665 * t844 + t666 * t842 + t667 * t840) * t916 + t931 * t679 + t930 * t678 + t932 * t677, (t665 * t835 + t666 * t831 + t667 * t833) * t852 + (-t665 * t843 - t666 * t841 - t667 * t839) * t916 + t921 * t679 + t923 * t678 + t922 * t677, 0; t814, 0, 0, t814 * t803, t620 * t858 + t622 * t857 + t624 * t856 + (t671 * t885 + t672 * t884 + t673 * t883) * t802, (t671 * t844 + t672 * t842 + t673 * t840 + t830 * t918 + t832 * t917 + t834 * t919) * t802 + (t827 - t822) * t754 + (-t828 - t821) * t753 + (-t829 - t823) * t752, (-t671 * t843 - t672 * t841 - t673 * t839 + t831 * t918 + t833 * t917 + t835 * t919) * t802 + (t824 - t818 / 0.2e1) * t754 + (t825 + t820 / 0.2e1) * t753 + (t826 + t819 / 0.2e1) * t752, 0;];
tau_reg  = t1;
