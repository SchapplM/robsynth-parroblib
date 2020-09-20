% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRR1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRR1G1P1A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:43
% EndTime: 2019-05-03 14:47:43
% DurationCPUTime: 0.36s
% Computational Cost: add. (144->65), mult. (284->129), div. (45->4), fcn. (295->14), ass. (0->64)
t825 = xP(3);
t811 = sin(t825);
t812 = cos(t825);
t826 = koppelP(3,2);
t829 = koppelP(3,1);
t791 = t811 * t829 + t812 * t826;
t816 = legFrame(3,3);
t805 = sin(t816);
t808 = cos(t816);
t835 = t811 * t826 - t812 * t829;
t853 = t791 * t808 + t805 * t835;
t827 = koppelP(2,2);
t830 = koppelP(2,1);
t792 = t811 * t830 + t812 * t827;
t817 = legFrame(2,3);
t806 = sin(t817);
t809 = cos(t817);
t834 = t811 * t827 - t812 * t830;
t852 = t792 * t809 + t806 * t834;
t828 = koppelP(1,2);
t831 = koppelP(1,1);
t793 = t811 * t831 + t812 * t828;
t818 = legFrame(1,3);
t807 = sin(t818);
t810 = cos(t818);
t833 = t811 * t828 - t812 * t831;
t851 = t793 * t810 + t807 * t833;
t821 = sin(qJ(2,1));
t815 = 0.1e1 / t821;
t850 = t851 * t815;
t819 = sin(qJ(2,3));
t813 = 0.1e1 / t819;
t849 = t853 * t813;
t820 = sin(qJ(2,2));
t814 = 0.1e1 / t820;
t848 = t852 * t814;
t797 = t805 * g(1) - t808 * g(2);
t844 = t797 * t813;
t798 = t806 * g(1) - t809 * g(2);
t843 = t798 * t814;
t799 = t807 * g(1) - t810 * g(2);
t842 = t799 * t815;
t841 = t805 * t813;
t840 = t806 * t814;
t839 = t807 * t815;
t838 = t808 * t813;
t837 = t809 * t814;
t836 = t810 * t815;
t832 = 0.1e1 / pkin(2);
t824 = cos(qJ(2,1));
t823 = cos(qJ(2,2));
t822 = cos(qJ(2,3));
t804 = t812 * g(1) + t811 * g(2);
t803 = t811 * g(1) - t812 * g(2);
t802 = t810 * g(1) + t807 * g(2);
t801 = t809 * g(1) + t806 * g(2);
t800 = t808 * g(1) + t805 * g(2);
t790 = -t799 * t821 + t802 * t824;
t789 = -t798 * t820 + t801 * t823;
t788 = -t797 * t819 + t800 * t822;
t787 = t799 * t824 + t802 * t821;
t786 = t798 * t823 + t801 * t820;
t785 = t797 * t822 + t800 * t819;
t1 = [((-t807 * t821 + t810 * t824) * t842 + (-t806 * t820 + t809 * t823) * t843 + (-t805 * t819 + t808 * t822) * t844) * MDP(1) + (-t811 * t803 - t812 * t804) * MDP(8) + ((-t785 * t838 - t786 * t837 - t787 * t836) * MDP(3) + (-t788 * t838 - t789 * t837 - t790 * t836) * MDP(4)) * t832; ((t807 * t824 + t821 * t810) * t842 + (t806 * t823 + t820 * t809) * t843 + (t805 * t822 + t819 * t808) * t844) * MDP(1) + (t812 * t803 - t811 * t804) * MDP(8) + ((-t785 * t841 - t786 * t840 - t787 * t839) * MDP(3) + (-t788 * t841 - t789 * t840 - t790 * t839) * MDP(4)) * t832; (((t807 * t793 - t810 * t833) * t821 - t851 * t824) * t842 + ((t806 * t792 - t809 * t834) * t820 - t852 * t823) * t843 + ((t805 * t791 - t808 * t835) * t819 - t853 * t822) * t844) * MDP(1) + t803 * MDP(6) + t804 * MDP(7) + ((t785 * t849 + t786 * t848 + t787 * t850) * MDP(3) + (t788 * t849 + t789 * t848 + t790 * t850) * MDP(4)) * t832;];
taugX  = t1;
