% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G1A0_invdyn_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:56
% EndTime: 2020-03-09 21:14:57
% DurationCPUTime: 1.71s
% Computational Cost: add. (16655->211), mult. (11100->369), div. (1944->9), fcn. (11400->30), ass. (0->194)
t786 = pkin(7) + qJ(2,1);
t776 = qJ(3,1) + t786;
t764 = sin(t776);
t767 = cos(t776);
t789 = legFrame(1,3);
t779 = sin(t789);
t782 = cos(t789);
t724 = t782 * t764 + t779 * t767;
t770 = sin(t786);
t773 = cos(t786);
t716 = t764 * t773 - t770 * t767;
t873 = 0.1e1 / t716;
t841 = t873 * t724;
t725 = -t764 * t779 + t782 * t767;
t840 = t873 * t725;
t785 = pkin(7) + qJ(2,2);
t775 = qJ(3,2) + t785;
t763 = sin(t775);
t766 = cos(t775);
t788 = legFrame(2,3);
t778 = sin(t788);
t781 = cos(t788);
t722 = t781 * t763 + t778 * t766;
t769 = sin(t785);
t772 = cos(t785);
t715 = t763 * t772 - t769 * t766;
t874 = 0.1e1 / t715;
t843 = t874 * t722;
t723 = -t763 * t778 + t781 * t766;
t842 = t874 * t723;
t784 = pkin(7) + qJ(2,3);
t774 = qJ(3,3) + t784;
t762 = sin(t774);
t765 = cos(t774);
t787 = legFrame(3,3);
t777 = sin(t787);
t780 = cos(t787);
t720 = t780 * t762 + t777 * t765;
t768 = sin(t784);
t771 = cos(t784);
t714 = t762 * t771 - t768 * t765;
t875 = 0.1e1 / t714;
t845 = t875 * t720;
t721 = -t762 * t777 + t780 * t765;
t844 = t875 * t721;
t877 = -2 * pkin(2);
t876 = 2 * pkin(2);
t802 = 1 / pkin(2);
t801 = 0.1e1 / pkin(3);
t798 = xDP(2);
t799 = xDP(1);
t726 = t777 * t798 + t780 * t799;
t727 = t777 * t799 - t780 * t798;
t681 = -t726 * t765 + t727 * t762;
t866 = (-(t726 * t771 - t768 * t727) * pkin(2) + t681 * pkin(3)) * t875;
t827 = t801 * t866;
t863 = t681 * t875;
t672 = (t827 - t863) * t802;
t806 = t762 * t768 + t765 * t771;
t666 = pkin(3) * t672 - t806 * t863;
t803 = 1 / pkin(2) ^ 2;
t830 = t672 * t866;
t791 = xDDP(1);
t838 = t791 * t802;
t790 = xDDP(2);
t839 = t790 * t802;
t837 = t838 * t844 + t839 * t845;
t709 = 0.1e1 / t714 ^ 2;
t862 = t681 * t709;
t660 = (-t666 * t862 + t830 * t875) * t803 + t837;
t872 = pkin(2) * t660;
t728 = t778 * t798 + t781 * t799;
t729 = t778 * t799 - t781 * t798;
t682 = -t728 * t766 + t729 * t763;
t865 = (-(t728 * t772 - t769 * t729) * pkin(2) + t682 * pkin(3)) * t874;
t826 = t801 * t865;
t861 = t682 * t874;
t673 = (t826 - t861) * t802;
t805 = t763 * t769 + t766 * t772;
t667 = pkin(3) * t673 - t805 * t861;
t829 = t673 * t865;
t836 = t838 * t842 + t839 * t843;
t711 = 0.1e1 / t715 ^ 2;
t860 = t682 * t711;
t661 = (-t667 * t860 + t829 * t874) * t803 + t836;
t871 = pkin(2) * t661;
t730 = t779 * t798 + t782 * t799;
t731 = t779 * t799 - t782 * t798;
t683 = -t730 * t767 + t731 * t764;
t864 = (-(t730 * t773 - t770 * t731) * pkin(2) + t683 * pkin(3)) * t873;
t825 = t801 * t864;
t859 = t683 * t873;
t674 = (t825 - t859) * t802;
t804 = t764 * t770 + t767 * t773;
t668 = t674 * pkin(3) - t804 * t859;
t828 = t674 * t864;
t835 = t838 * t840 + t839 * t841;
t713 = 0.1e1 / t716 ^ 2;
t858 = t683 * t713;
t662 = (-t668 * t858 + t828 * t873) * t803 + t835;
t870 = pkin(2) * t662;
t669 = (-t863 + t827 / 0.2e1) * t802;
t800 = pkin(3) ^ 2;
t831 = 0.2e1 * pkin(3);
t869 = (-t672 * t800 + (-t806 * t669 * t831 + t863) * pkin(2)) * t801;
t670 = (-t861 + t826 / 0.2e1) * t802;
t868 = (-t673 * t800 + (-t805 * t670 * t831 + t861) * pkin(2)) * t801;
t671 = (-t859 + t825 / 0.2e1) * t802;
t867 = (-t674 * t800 + (-t804 * t671 * t831 + t859) * pkin(2)) * t801;
t690 = pkin(2) * (t768 * t780 + t777 * t771) + t720 * pkin(3);
t857 = t690 * t875;
t856 = t690 * t790;
t691 = pkin(2) * (t769 * t781 + t778 * t772) + t722 * pkin(3);
t855 = t691 * t874;
t854 = t691 * t790;
t692 = pkin(2) * (t770 * t782 + t779 * t773) + t724 * pkin(3);
t853 = t692 * t873;
t852 = t692 * t790;
t693 = -pkin(2) * (t768 * t777 - t780 * t771) + t721 * pkin(3);
t851 = t693 * t875;
t850 = t693 * t791;
t694 = -pkin(2) * (t769 * t778 - t781 * t772) + t723 * pkin(3);
t849 = t694 * t874;
t848 = t694 * t791;
t695 = -pkin(2) * (t770 * t779 - t782 * t773) + t725 * pkin(3);
t847 = t695 * t873;
t846 = t695 * t791;
t759 = t787 + t774;
t753 = sin(t759);
t756 = cos(t759);
t834 = g(1) * t756 + g(2) * t753;
t760 = t788 + t775;
t754 = sin(t760);
t757 = cos(t760);
t833 = g(1) * t757 + g(2) * t754;
t761 = t789 + t776;
t755 = sin(t761);
t758 = cos(t761);
t832 = g(1) * t758 + g(2) * t755;
t824 = t681 ^ 2 * t709 * t802;
t823 = t682 ^ 2 * t711 * t802;
t822 = t683 ^ 2 * t713 * t802;
t821 = t803 * t862;
t820 = t803 * t860;
t819 = t803 * t858;
t818 = g(1) * t753 - g(2) * t756;
t817 = g(1) * t754 - g(2) * t757;
t816 = g(1) * t755 - g(2) * t758;
t815 = t803 * t830;
t814 = t803 * t829;
t813 = t803 * t828;
t812 = t669 * t802 * t827;
t811 = t670 * t802 * t826;
t810 = t671 * t802 * t825;
t809 = (t806 * pkin(2) + pkin(3)) * t815;
t808 = (t805 * pkin(2) + pkin(3)) * t814;
t807 = (t804 * pkin(2) + pkin(3)) * t813;
t797 = cos(qJ(3,1));
t796 = cos(qJ(3,2));
t795 = cos(qJ(3,3));
t794 = sin(qJ(3,1));
t793 = sin(qJ(3,2));
t792 = sin(qJ(3,3));
t737 = t782 * g(1) + t779 * g(2);
t736 = t781 * g(1) + t778 * g(2);
t735 = t780 * g(1) + t777 * g(2);
t734 = t779 * g(1) - t782 * g(2);
t733 = t778 * g(1) - t781 * g(2);
t732 = t777 * g(1) - t780 * g(2);
t701 = -t734 * t770 + t737 * t773;
t700 = -t733 * t769 + t736 * t772;
t699 = -t732 * t768 + t735 * t771;
t698 = t734 * t773 + t737 * t770;
t697 = t733 * t772 + t736 * t769;
t696 = t732 * t771 + t735 * t768;
t659 = t793 * t823 + t796 * t871 + t817;
t658 = t792 * t824 + t795 * t872 + t818;
t657 = -t794 * t870 + t797 * t822 + t832;
t656 = -t793 * t871 + t796 * t823 + t833;
t655 = -t792 * t872 + t795 * t824 + t834;
t654 = t794 * t822 + t797 * t870 + t816;
t653 = (-t668 - t867) * t819 - (-t813 + (t807 + (t846 + t852) * t802) * t801) * t873 + t835;
t652 = (-t667 - t868) * t820 - (-t814 + (t808 + (t848 + t854) * t802) * t801) * t874 + t836;
t651 = (-t666 - t869) * t821 - (-t815 + (t809 + (t850 + t856) * t802) * t801) * t875 + t837;
t650 = (-t668 - t867 / 0.2e1) * t819 - (-t813 + (t807 / 0.2e1 + (t846 / 0.2e1 + t852 / 0.2e1) * t802) * t801) * t873 + t835;
t649 = (-t667 - t868 / 0.2e1) * t820 - (-t814 + (t808 / 0.2e1 + (t848 / 0.2e1 + t854 / 0.2e1) * t802) * t801) * t874 + t836;
t648 = (-t666 - t869 / 0.2e1) * t821 - (-t815 + (t809 / 0.2e1 + (t850 / 0.2e1 + t856 / 0.2e1) * t802) * t801) * t875 + t837;
t647 = (t794 * t650 + t797 * t810) * t877 + t832;
t646 = (t793 * t649 + t796 * t811) * t877 + t833;
t645 = (t792 * t648 + t795 * t812) * t877 + t834;
t644 = (t650 * t797 - t794 * t810) * t876 + t816;
t643 = (t649 * t796 - t793 * t811) * t876 + t817;
t642 = (t648 * t795 - t792 * t812) * t876 + t818;
t1 = [(t791 - g(1)) * MDP(8) + ((t660 * t844 + t661 * t842 + t662 * t840) * MDP(2) + (t696 * t844 + t697 * t842 + t698 * t840) * MDP(3) + (t699 * t844 + t700 * t842 + t701 * t840) * MDP(4) + (t651 * t844 + t652 * t842 + t653 * t840) * MDP(5) + (t642 * t844 + t643 * t842 + t644 * t840) * MDP(6) + (t645 * t844 + t646 * t842 + t647 * t840) * MDP(7) + ((-t651 * t851 - t652 * t849 - t653 * t847) * MDP(5) + (-t654 * t847 - t658 * t851 - t659 * t849) * MDP(6) + (-t655 * t851 - t656 * t849 - t657 * t847) * MDP(7)) * t801) * t802; (t790 - g(2)) * MDP(8) + ((t660 * t845 + t661 * t843 + t662 * t841) * MDP(2) + (t696 * t845 + t697 * t843 + t698 * t841) * MDP(3) + (t699 * t845 + t700 * t843 + t701 * t841) * MDP(4) + (t651 * t845 + t652 * t843 + t653 * t841) * MDP(5) + (t642 * t845 + t643 * t843 + t644 * t841) * MDP(6) + (t645 * t845 + t646 * t843 + t647 * t841) * MDP(7) + ((-t651 * t857 - t652 * t855 - t653 * t853) * MDP(5) + (-t654 * t853 - t658 * t857 - t659 * t855) * MDP(6) + (-t655 * t857 - t656 * t855 - t657 * t853) * MDP(7)) * t801) * t802; ((3 * MDP(1)) + MDP(8)) * (xDDP(3) - g(3));];
tauX  = t1;
