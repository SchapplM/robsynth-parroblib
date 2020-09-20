% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR2G3P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%   see P3PRRR2G3P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR2G3P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR2G3P1A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:10
% EndTime: 2020-03-09 21:20:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (289->70), mult. (289->116), div. (84->5), fcn. (300->33), ass. (0->72)
t1029 = legFrame(3,2);
t1054 = sin(t1029);
t1030 = legFrame(2,2);
t1053 = sin(t1030);
t1031 = legFrame(1,2);
t1052 = sin(t1031);
t1011 = -t1029 + qJ(2,3);
t1023 = 0.1e1 / sin(qJ(3,3));
t1005 = qJ(3,3) + t1011;
t999 = sin(t1005);
t1051 = t1023 * (pkin(2) * t999 + pkin(1) * sin(t1011));
t1002 = cos(t1005);
t1050 = t1023 * (-pkin(2) * t1002 - pkin(1) * cos(t1011));
t1049 = t1023 * t999;
t1012 = -t1030 + qJ(2,2);
t1006 = qJ(3,2) + t1012;
t1000 = sin(t1006);
t1024 = 0.1e1 / sin(qJ(3,2));
t1048 = t1024 * (pkin(2) * t1000 + pkin(1) * sin(t1012));
t1003 = cos(t1006);
t1047 = t1024 * (-pkin(2) * t1003 - pkin(1) * cos(t1012));
t1013 = -t1031 + qJ(2,1);
t1007 = qJ(3,1) + t1013;
t1001 = sin(t1007);
t1025 = 0.1e1 / sin(qJ(3,1));
t1046 = t1025 * (pkin(2) * t1001 + pkin(1) * sin(t1013));
t1004 = cos(t1007);
t1045 = t1025 * (-pkin(2) * t1004 - pkin(1) * cos(t1013));
t1044 = t1000 * t1024;
t1043 = t1001 * t1025;
t1042 = t1002 * t1023;
t1041 = t1003 * t1024;
t1040 = t1004 * t1025;
t1039 = 0.1e1 / pkin(1);
t1038 = 0.1e1 / pkin(2);
t1037 = cos(qJ(2,1));
t1036 = cos(qJ(2,2));
t1035 = cos(qJ(2,3));
t1034 = sin(qJ(2,1));
t1033 = sin(qJ(2,2));
t1032 = sin(qJ(2,3));
t1028 = qJ(2,1) + qJ(3,1);
t1027 = qJ(2,2) + qJ(3,2);
t1026 = qJ(2,3) + qJ(3,3);
t1022 = cos(t1031);
t1021 = cos(t1030);
t1020 = cos(t1029);
t1019 = cos(t1028);
t1018 = cos(t1027);
t1017 = cos(t1026);
t1016 = sin(t1028);
t1015 = sin(t1027);
t1014 = sin(t1026);
t998 = t1022 * g(1) - t1052 * g(2);
t997 = t1021 * g(1) - t1053 * g(2);
t996 = t1020 * g(1) - t1054 * g(2);
t995 = t1052 * g(1) + t1022 * g(2);
t994 = t1053 * g(1) + t1021 * g(2);
t993 = t1054 * g(1) + t1020 * g(2);
t986 = -t998 * t1034 + t995 * t1037;
t985 = t995 * t1034 + t998 * t1037;
t984 = -t997 * t1033 + t994 * t1036;
t983 = t994 * t1033 + t997 * t1036;
t982 = -t996 * t1032 + t993 * t1035;
t981 = t993 * t1032 + t996 * t1035;
t980 = -t998 * t1016 + t995 * t1019;
t979 = t995 * t1016 + t998 * t1019;
t978 = -t997 * t1015 + t994 * t1018;
t977 = t994 * t1015 + t997 * t1018;
t976 = -t996 * t1014 + t993 * t1017;
t975 = t993 * t1014 + t996 * t1017;
t1 = [-g(1) * MDP(8) + ((-t985 * t1043 - t983 * t1044 - t981 * t1049) * MDP(3) + (-t986 * t1043 - t984 * t1044 - t982 * t1049) * MDP(4) + (-t979 * t1043 - t977 * t1044 - t975 * t1049) * MDP(6) + (-t980 * t1043 - t978 * t1044 - t976 * t1049) * MDP(7) + ((t979 * t1046 + t977 * t1048 + t975 * t1051) * MDP(6) + (t980 * t1046 + t978 * t1048 + t976 * t1051) * MDP(7)) * t1038) * t1039; -g(2) * MDP(8) + ((t985 * t1040 + t983 * t1041 + t981 * t1042) * MDP(3) + (t986 * t1040 + t984 * t1041 + t982 * t1042) * MDP(4) + (t979 * t1040 + t977 * t1041 + t975 * t1042) * MDP(6) + (t980 * t1040 + t978 * t1041 + t976 * t1042) * MDP(7) + ((t979 * t1045 + t977 * t1047 + t975 * t1050) * MDP(6) + (t980 * t1045 + t978 * t1047 + t976 * t1050) * MDP(7)) * t1038) * t1039; (-(3 * MDP(1)) - MDP(8)) * g(3);];
taugX  = t1;
