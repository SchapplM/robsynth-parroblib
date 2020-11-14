% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR2G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR2G3A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRRR2G3A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G3A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:08
% EndTime: 2020-03-09 21:20:09
% DurationCPUTime: 0.51s
% Computational Cost: add. (719->105), mult. (353->197), div. (408->9), fcn. (486->18), ass. (0->96)
t822 = MDP(2) / pkin(1) ^ 2;
t758 = -legFrame(3,2) + qJ(2,3);
t755 = qJ(3,3) + t758;
t749 = sin(t755);
t743 = pkin(2) * t749 + pkin(1) * sin(t758);
t770 = sin(qJ(3,3));
t761 = 0.1e1 / t770;
t821 = t743 * t761;
t759 = -legFrame(2,2) + qJ(2,2);
t756 = qJ(3,2) + t759;
t750 = sin(t756);
t744 = pkin(2) * t750 + pkin(1) * sin(t759);
t771 = sin(qJ(3,2));
t763 = 0.1e1 / t771;
t820 = t744 * t763;
t760 = -legFrame(1,2) + qJ(2,1);
t757 = qJ(3,1) + t760;
t751 = sin(t757);
t745 = pkin(2) * t751 + pkin(1) * sin(t760);
t772 = sin(qJ(3,1));
t765 = 0.1e1 / t772;
t819 = t745 * t765;
t752 = cos(t755);
t746 = -pkin(2) * t752 - pkin(1) * cos(t758);
t818 = t746 * t761;
t753 = cos(t756);
t747 = -pkin(2) * t753 - pkin(1) * cos(t759);
t817 = t747 * t763;
t754 = cos(t757);
t748 = -pkin(2) * t754 - pkin(1) * cos(t760);
t816 = t748 * t765;
t815 = t749 * t761;
t814 = t750 * t763;
t813 = t751 * t765;
t812 = t752 * t761;
t811 = t753 * t763;
t810 = t754 * t765;
t773 = cos(qJ(3,3));
t809 = t761 * t773;
t777 = 0.1e1 / pkin(1);
t808 = t761 * t777;
t762 = 0.1e1 / t770 ^ 2;
t807 = t762 * t773;
t774 = cos(qJ(3,2));
t806 = t763 * t774;
t805 = t763 * t777;
t764 = 0.1e1 / t771 ^ 2;
t804 = t764 * t774;
t775 = cos(qJ(3,1));
t803 = t765 * t775;
t802 = t765 * t777;
t766 = 0.1e1 / t772 ^ 2;
t801 = t766 * t775;
t776 = 0.1e1 / pkin(2);
t800 = t776 * t777;
t799 = t749 * t809;
t798 = t749 * t807;
t797 = t750 * t806;
t796 = t750 * t804;
t795 = t751 * t803;
t794 = t751 * t801;
t793 = t752 * t809;
t792 = t752 * t807;
t791 = t753 * t806;
t790 = t753 * t804;
t789 = t754 * t803;
t788 = t754 * t801;
t787 = t761 * t800;
t786 = t763 * t800;
t785 = t765 * t800;
t784 = t749 * t808;
t783 = t750 * t805;
t782 = t751 * t802;
t781 = t752 * t808;
t780 = t753 * t805;
t779 = t754 * t802;
t742 = t748 * t785;
t741 = t747 * t786;
t740 = t746 * t787;
t739 = t745 * t785;
t738 = t744 * t786;
t737 = t743 * t787;
t736 = t742 + t779;
t735 = t741 + t780;
t734 = t740 + t781;
t733 = t739 - t782;
t732 = t738 - t783;
t731 = t737 - t784;
t730 = t742 + 0.2e1 * t779;
t729 = t741 + 0.2e1 * t780;
t728 = t740 + 0.2e1 * t781;
t727 = t739 - 0.2e1 * t782;
t726 = t738 - 0.2e1 * t783;
t725 = t737 - 0.2e1 * t784;
t724 = (-t749 * t752 * t762 - t750 * t753 * t764 - t751 * t754 * t766) * t822;
t1 = [(-t725 * t799 - t726 * t797 - t727 * t795) * MDP(6) + (t749 * t725 + t750 * t726 + t751 * t727) * MDP(7) + MDP(8) + (t749 ^ 2 * t762 + t750 ^ 2 * t764 + t751 ^ 2 * t766) * t822 + ((-t731 * t815 - t732 * t814 - t733 * t813) * MDP(5) + ((t731 * t821 + t732 * t820 + t733 * t819) * MDP(5) + (-t743 * t798 - t744 * t796 - t745 * t794) * MDP(6) + (t743 * t815 + t744 * t814 + t745 * t813) * MDP(7)) * t776) * t777; t724 + (-t728 * t799 - t729 * t797 - t730 * t795) * MDP(6) + (t749 * t728 + t750 * t729 + t751 * t730) * MDP(7) + ((-t734 * t815 - t735 * t814 - t736 * t813) * MDP(5) + ((t734 * t821 + t735 * t820 + t736 * t819) * MDP(5) + (t743 * t792 + t744 * t790 + t745 * t788) * MDP(6) + (-t743 * t812 - t744 * t811 - t745 * t810) * MDP(7)) * t776) * t777; 0; t724 + (t725 * t793 + t726 * t791 + t727 * t789) * MDP(6) + (-t752 * t725 - t753 * t726 - t754 * t727) * MDP(7) + ((t731 * t812 + t732 * t811 + t733 * t810) * MDP(5) + ((t731 * t818 + t732 * t817 + t733 * t816) * MDP(5) + (-t746 * t798 - t747 * t796 - t748 * t794) * MDP(6) + (t746 * t815 + t747 * t814 + t748 * t813) * MDP(7)) * t776) * t777; (t728 * t793 + t729 * t791 + t730 * t789) * MDP(6) + (-t752 * t728 - t753 * t729 - t754 * t730) * MDP(7) + MDP(8) + (t752 ^ 2 * t762 + t753 ^ 2 * t764 + t754 ^ 2 * t766) * t822 + ((t734 * t812 + t735 * t811 + t736 * t810) * MDP(5) + ((t734 * t818 + t735 * t817 + t736 * t816) * MDP(5) + (t746 * t792 + t747 * t790 + t748 * t788) * MDP(6) + (-t746 * t812 - t747 * t811 - t748 * t810) * MDP(7)) * t776) * t777; 0; 0; 0; (3 * MDP(1)) + MDP(8);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
