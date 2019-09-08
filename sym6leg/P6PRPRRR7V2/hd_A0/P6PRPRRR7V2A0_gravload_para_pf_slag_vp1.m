% Calculate Gravitation load for parallel robot
% P6PRPRRR7V2G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d2,d4,theta1,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-17 04:21
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(10,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PRPRRR7V2G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-17 02:52:42
% EndTime: 2019-05-17 02:53:02
% DurationCPUTime: 22.47s
% Computational Cost: add. (11813->738), mult. (26241->1382), div. (162->18), fcn. (26501->64), ass. (0->561)
t5279 = cos(pkin(6)) * pkin(8);
t4979 = qJ(3,3) + t5279;
t5034 = cos(pkin(5));
t5025 = t5034 ^ 2;
t5035 = cos(pkin(4));
t5030 = sin(pkin(4));
t5162 = t5030 * t5034;
t5031 = cos(pkin(10));
t5272 = sin(pkin(10));
t5280 = pkin(8) * sin(pkin(6));
t4954 = -pkin(3) * t5272 + t5031 * t5280;
t4935 = t4954 * t5034;
t5029 = sin(pkin(5));
t5309 = -t5029 * t4935 - t5279;
t5190 = t4979 * t5029;
t4879 = -t4935 - t5190;
t4938 = pkin(3) * t5031 + t5272 * t5280 + pkin(2);
t5060 = sin(qJ(2,3));
t4916 = t5060 * t4938;
t5063 = cos(qJ(2,3));
t5346 = t4879 * t5063 + t4916;
t4705 = 0.1e1 / (t5162 * t5346 + (t4979 * t5025 - qJ(3,3) + t5309) * t5035);
t4980 = qJ(3,2) + t5279;
t5189 = t4980 * t5029;
t4880 = -t4935 - t5189;
t5061 = sin(qJ(2,2));
t4917 = t5061 * t4938;
t5064 = cos(qJ(2,2));
t5345 = t4880 * t5064 + t4917;
t4706 = 0.1e1 / (t5162 * t5345 + (t4980 * t5025 - qJ(3,2) + t5309) * t5035);
t4981 = qJ(3,1) + t5279;
t5188 = t4981 * t5029;
t4881 = -t4935 - t5188;
t5062 = sin(qJ(2,1));
t4918 = t5062 * t4938;
t5065 = cos(qJ(2,1));
t5344 = t4881 * t5065 + t4918;
t4707 = 0.1e1 / (t5162 * t5344 + (t4981 * t5025 - qJ(3,1) + t5309) * t5035);
t4976 = qJ(3,6) + t5279;
t5193 = t4976 * t5029;
t4876 = -t4935 - t5193;
t5048 = sin(qJ(2,6));
t4909 = t5048 * t4938;
t5051 = cos(qJ(2,6));
t5349 = t4876 * t5051 + t4909;
t4702 = 0.1e1 / (t5162 * t5349 + (t4976 * t5025 - qJ(3,6) + t5309) * t5035);
t4977 = qJ(3,5) + t5279;
t5192 = t4977 * t5029;
t4877 = -t4935 - t5192;
t5049 = sin(qJ(2,5));
t4910 = t5049 * t4938;
t5052 = cos(qJ(2,5));
t5348 = t4877 * t5052 + t4910;
t4703 = 0.1e1 / (t5162 * t5348 + (t4977 * t5025 - qJ(3,5) + t5309) * t5035);
t4978 = qJ(3,4) + t5279;
t5191 = t4978 * t5029;
t4878 = -t4935 - t5191;
t5050 = sin(qJ(2,4));
t4911 = t5050 * t4938;
t5053 = cos(qJ(2,4));
t5347 = t4878 * t5053 + t4911;
t4704 = 0.1e1 / (t5162 * t5347 + (t4978 * t5025 - qJ(3,4) + t5309) * t5035);
t4860 = -t5065 * t5188 + t4918;
t4934 = t5035 * t4954;
t4908 = t4934 * t5029;
t5138 = t4860 * t5030 - t4908;
t4972 = t5035 * t4981;
t5156 = t5030 * t5065;
t5305 = t4954 * t5156 - t4972;
t5144 = -t5025 * t5305 - t4972;
t5343 = 0.1e1 / (t5034 * t5138 + t5144);
t4859 = -t5064 * t5189 + t4917;
t5139 = t4859 * t5030 - t4908;
t4971 = t5035 * t4980;
t5157 = t5030 * t5064;
t5304 = t4954 * t5157 - t4971;
t5145 = -t5025 * t5304 - t4971;
t5342 = 0.1e1 / (t5034 * t5139 + t5145);
t4858 = -t5063 * t5190 + t4916;
t5140 = t4858 * t5030 - t4908;
t4970 = t5035 * t4979;
t5158 = t5030 * t5063;
t5303 = t4954 * t5158 - t4970;
t5146 = -t5025 * t5303 - t4970;
t5341 = 0.1e1 / (t5034 * t5140 + t5146);
t4854 = -t5053 * t5191 + t4911;
t5141 = t4854 * t5030 - t4908;
t4966 = t5035 * t4978;
t5159 = t5030 * t5053;
t5302 = t4954 * t5159 - t4966;
t5147 = -t5025 * t5302 - t4966;
t5340 = 0.1e1 / (t5034 * t5141 + t5147);
t4853 = -t5052 * t5192 + t4910;
t5142 = t4853 * t5030 - t4908;
t4965 = t5035 * t4977;
t5160 = t5030 * t5052;
t5301 = t4954 * t5160 - t4965;
t5148 = -t5025 * t5301 - t4965;
t5339 = 0.1e1 / (t5034 * t5142 + t5148);
t4852 = -t5051 * t5193 + t4909;
t5143 = t4852 * t5030 - t4908;
t4964 = t5035 * t4976;
t5161 = t5030 * t5051;
t5300 = t4954 * t5161 - t4964;
t5149 = -t5025 * t5300 - t4964;
t5338 = 0.1e1 / (t5034 * t5143 + t5149);
t5068 = xP(5);
t5020 = sin(t5068);
t5023 = cos(t5068);
t5073 = koppelP(6,3);
t5067 = xP(6);
t5019 = sin(t5067);
t5022 = cos(t5067);
t5079 = koppelP(6,1);
t5288 = koppelP(6,2);
t5299 = t5019 * t5288 - t5022 * t5079;
t4894 = t5020 * t5299 + t5023 * t5073;
t4945 = t5019 * t5079 + t5022 * t5288;
t5069 = xP(4);
t5021 = sin(t5069);
t5024 = cos(t5069);
t4780 = t4894 * t5021 - t4945 * t5024;
t4781 = t4894 * t5024 + t4945 * t5021;
t5074 = koppelP(5,3);
t5080 = koppelP(5,1);
t5289 = koppelP(5,2);
t5298 = t5019 * t5289 - t5022 * t5080;
t4896 = t5020 * t5298 + t5023 * t5074;
t4946 = t5019 * t5080 + t5022 * t5289;
t4784 = t4896 * t5021 - t4946 * t5024;
t4785 = t4896 * t5024 + t4946 * t5021;
t5075 = koppelP(4,3);
t5081 = koppelP(4,1);
t5290 = koppelP(4,2);
t5297 = t5019 * t5290 - t5022 * t5081;
t4898 = t5020 * t5297 + t5023 * t5075;
t4947 = t5019 * t5081 + t5022 * t5290;
t4788 = t4898 * t5021 - t4947 * t5024;
t4789 = t4898 * t5024 + t4947 * t5021;
t5076 = koppelP(3,3);
t5082 = koppelP(3,1);
t5291 = koppelP(3,2);
t5296 = t5019 * t5291 - t5022 * t5082;
t4900 = t5020 * t5296 + t5023 * t5076;
t4948 = t5019 * t5082 + t5022 * t5291;
t4792 = t4900 * t5021 - t4948 * t5024;
t4793 = t4900 * t5024 + t4948 * t5021;
t5077 = koppelP(2,3);
t5083 = koppelP(2,1);
t5292 = koppelP(2,2);
t5295 = t5019 * t5292 - t5022 * t5083;
t4902 = t5020 * t5295 + t5023 * t5077;
t4949 = t5019 * t5083 + t5022 * t5292;
t4796 = t4902 * t5021 - t4949 * t5024;
t4797 = t4902 * t5024 + t4949 * t5021;
t5078 = koppelP(1,3);
t5084 = koppelP(1,1);
t5293 = koppelP(1,2);
t5294 = t5019 * t5293 - t5022 * t5084;
t4904 = t5020 * t5294 + t5023 * t5078;
t4950 = t5019 * t5084 + t5022 * t5293;
t4800 = t4904 * t5021 - t4950 * t5024;
t4801 = t4904 * t5024 + t4950 * t5021;
t5070 = rSges(4,3);
t5071 = rSges(4,2);
t5072 = rSges(4,1);
t5098 = t5019 * t5071 - t5022 * t5072;
t5331 = -t5020 * t5070 + t5023 * t5098;
t5036 = legFrame(6,3);
t4983 = sin(t5036);
t4995 = cos(t5036);
t5027 = sin(pkin(9));
t5032 = cos(pkin(9));
t4922 = t4983 * t5032 + t4995 * t5027;
t5321 = t4922 * t5030;
t5037 = legFrame(5,3);
t4984 = sin(t5037);
t4996 = cos(t5037);
t4923 = t4984 * t5032 + t4996 * t5027;
t5320 = t4923 * t5030;
t5038 = legFrame(4,3);
t4985 = sin(t5038);
t4997 = cos(t5038);
t4924 = t4985 * t5032 + t4997 * t5027;
t5319 = t4924 * t5030;
t5039 = legFrame(3,3);
t4986 = sin(t5039);
t4998 = cos(t5039);
t4925 = t4986 * t5032 + t4998 * t5027;
t5318 = t4925 * t5030;
t5040 = legFrame(2,3);
t4987 = sin(t5040);
t4999 = cos(t5040);
t4926 = t4987 * t5032 + t4999 * t5027;
t5317 = t4926 * t5030;
t5041 = legFrame(1,3);
t4988 = sin(t5041);
t5000 = cos(t5041);
t4927 = t4988 * t5032 + t5000 * t5027;
t5316 = t4927 * t5030;
t5042 = legFrame(6,1);
t4989 = sin(t5042);
t5001 = cos(t5042);
t5054 = legFrame(6,2);
t5007 = sin(t5054);
t5181 = t5001 * t5007;
t5187 = t4989 * t5007;
t5013 = cos(t5054);
t5278 = g(1) * t5013;
t4732 = -t4983 * t5278 + (-t4983 * t5187 + t4995 * t5001) * g(2) + (t4983 * t5181 + t4989 * t4995) * g(3);
t4733 = t4995 * t5278 + (t4983 * t5001 + t4995 * t5187) * g(2) + (t4983 * t4989 - t4995 * t5181) * g(3);
t4666 = -t5032 * t4732 + t4733 * t5027;
t4882 = t5007 * g(1) + (-g(2) * t4989 + g(3) * t5001) * t5013;
t5104 = t4732 * t5027 + t4733 * t5032;
t5110 = t4666 * t5035 - t4882 * t5030;
t5090 = t5048 * t5104 + t5051 * t5110;
t4606 = t5090 * t5029 + (t4666 * t5030 + t4882 * t5035) * t5034;
t5287 = m(3) * t4606;
t5043 = legFrame(5,1);
t4990 = sin(t5043);
t5002 = cos(t5043);
t5055 = legFrame(5,2);
t5008 = sin(t5055);
t5180 = t5002 * t5008;
t5186 = t4990 * t5008;
t5014 = cos(t5055);
t5277 = g(1) * t5014;
t4734 = -t4984 * t5277 + (-t4984 * t5186 + t4996 * t5002) * g(2) + (t4984 * t5180 + t4990 * t4996) * g(3);
t4735 = t4996 * t5277 + (t4984 * t5002 + t4996 * t5186) * g(2) + (t4984 * t4990 - t4996 * t5180) * g(3);
t4667 = -t5032 * t4734 + t4735 * t5027;
t4883 = t5008 * g(1) + (-g(2) * t4990 + g(3) * t5002) * t5014;
t5103 = t4734 * t5027 + t4735 * t5032;
t5109 = t4667 * t5035 - t4883 * t5030;
t5089 = t5049 * t5103 + t5052 * t5109;
t4607 = t5089 * t5029 + (t4667 * t5030 + t4883 * t5035) * t5034;
t5286 = m(3) * t4607;
t5044 = legFrame(4,1);
t4991 = sin(t5044);
t5003 = cos(t5044);
t5056 = legFrame(4,2);
t5009 = sin(t5056);
t5179 = t5003 * t5009;
t5185 = t4991 * t5009;
t5015 = cos(t5056);
t5276 = g(1) * t5015;
t4736 = -t4985 * t5276 + (-t4985 * t5185 + t4997 * t5003) * g(2) + (t4985 * t5179 + t4991 * t4997) * g(3);
t4737 = t4997 * t5276 + (t4985 * t5003 + t4997 * t5185) * g(2) + (t4985 * t4991 - t4997 * t5179) * g(3);
t4668 = -t5032 * t4736 + t4737 * t5027;
t4884 = t5009 * g(1) + (-g(2) * t4991 + g(3) * t5003) * t5015;
t5102 = t4736 * t5027 + t4737 * t5032;
t5108 = t4668 * t5035 - t4884 * t5030;
t5088 = t5050 * t5102 + t5053 * t5108;
t4608 = t5088 * t5029 + (t4668 * t5030 + t4884 * t5035) * t5034;
t5285 = m(3) * t4608;
t5045 = legFrame(3,1);
t4992 = sin(t5045);
t5004 = cos(t5045);
t5057 = legFrame(3,2);
t5010 = sin(t5057);
t5178 = t5004 * t5010;
t5184 = t4992 * t5010;
t5016 = cos(t5057);
t5275 = g(1) * t5016;
t4738 = -t4986 * t5275 + (-t4986 * t5184 + t4998 * t5004) * g(2) + (t4986 * t5178 + t4992 * t4998) * g(3);
t4739 = t4998 * t5275 + (t4986 * t5004 + t4998 * t5184) * g(2) + (t4986 * t4992 - t4998 * t5178) * g(3);
t4669 = -t5032 * t4738 + t4739 * t5027;
t4885 = t5010 * g(1) + (-g(2) * t4992 + g(3) * t5004) * t5016;
t5101 = t4738 * t5027 + t4739 * t5032;
t5107 = t4669 * t5035 - t4885 * t5030;
t5087 = t5060 * t5101 + t5063 * t5107;
t4609 = t5087 * t5029 + (t4669 * t5030 + t4885 * t5035) * t5034;
t5284 = m(3) * t4609;
t5046 = legFrame(2,1);
t4993 = sin(t5046);
t5005 = cos(t5046);
t5058 = legFrame(2,2);
t5011 = sin(t5058);
t5177 = t5005 * t5011;
t5183 = t4993 * t5011;
t5017 = cos(t5058);
t5274 = g(1) * t5017;
t4740 = -t4987 * t5274 + (-t4987 * t5183 + t4999 * t5005) * g(2) + (t4987 * t5177 + t4993 * t4999) * g(3);
t4741 = t4999 * t5274 + (t4987 * t5005 + t4999 * t5183) * g(2) + (t4987 * t4993 - t4999 * t5177) * g(3);
t4670 = -t5032 * t4740 + t4741 * t5027;
t4886 = t5011 * g(1) + (-g(2) * t4993 + g(3) * t5005) * t5017;
t5100 = t4740 * t5027 + t4741 * t5032;
t5106 = t4670 * t5035 - t4886 * t5030;
t5086 = t5061 * t5100 + t5064 * t5106;
t4610 = t5086 * t5029 + (t4670 * t5030 + t4886 * t5035) * t5034;
t5283 = m(3) * t4610;
t5047 = legFrame(1,1);
t4994 = sin(t5047);
t5006 = cos(t5047);
t5059 = legFrame(1,2);
t5012 = sin(t5059);
t5176 = t5006 * t5012;
t5182 = t4994 * t5012;
t5018 = cos(t5059);
t5273 = g(1) * t5018;
t4742 = -t4988 * t5273 + (-t4988 * t5182 + t5000 * t5006) * g(2) + (t4988 * t5176 + t4994 * t5000) * g(3);
t4743 = t5000 * t5273 + (t4988 * t5006 + t5000 * t5182) * g(2) + (t4988 * t4994 - t5000 * t5176) * g(3);
t4671 = -t5032 * t4742 + t4743 * t5027;
t4887 = t5012 * g(1) + (-g(2) * t4994 + g(3) * t5006) * t5018;
t5099 = t4742 * t5027 + t4743 * t5032;
t5105 = t4671 * t5035 - t4887 * t5030;
t5085 = t5062 * t5099 + t5065 * t5105;
t4611 = t5085 * t5029 + (t4671 * t5030 + t4887 * t5035) * t5034;
t5282 = m(3) * t4611;
t5281 = m(3) * t5029;
t4915 = m(2) * rSges(2,1) + (rSges(3,1) * t5031 - rSges(3,2) * t5272 + pkin(2)) * m(3);
t5137 = m(3) * (rSges(3,1) * t5272 + t5031 * rSges(3,2)) * t5034 + m(2) * rSges(2,2);
t4600 = t5090 * t4915 + (-t5048 * t5110 + t5051 * t5104) * (-(qJ(3,6) + rSges(3,3)) * t5281 + t5137);
t4680 = 0.1e1 / (((-t4976 * t5161 - t4934) * t5029 + t5030 * t4909) * t5034 + t5149);
t5271 = t4600 * t4680;
t4601 = t5089 * t4915 + (-t5049 * t5109 + t5052 * t5103) * (-(qJ(3,5) + rSges(3,3)) * t5281 + t5137);
t4683 = 0.1e1 / (((-t4977 * t5160 - t4934) * t5029 + t5030 * t4910) * t5034 + t5148);
t5270 = t4601 * t4683;
t4602 = t5088 * t4915 + (-t5050 * t5108 + t5053 * t5102) * (-(qJ(3,4) + rSges(3,3)) * t5281 + t5137);
t4686 = 0.1e1 / (((-t4978 * t5159 - t4934) * t5029 + t5030 * t4911) * t5034 + t5147);
t5269 = t4602 * t4686;
t4603 = t5087 * t4915 + (-t5060 * t5107 + t5063 * t5101) * (-(qJ(3,3) + rSges(3,3)) * t5281 + t5137);
t4689 = 0.1e1 / (((-t4979 * t5158 - t4934) * t5029 + t5030 * t4916) * t5034 + t5146);
t5268 = t4603 * t4689;
t4604 = t5086 * t4915 + (-t5061 * t5106 + t5064 * t5100) * (-(qJ(3,2) + rSges(3,3)) * t5281 + t5137);
t4692 = 0.1e1 / (((-t4980 * t5157 - t4934) * t5029 + t5030 * t4917) * t5034 + t5145);
t5267 = t4604 * t4692;
t4605 = t5085 * t4915 + (-t5062 * t5105 + t5065 * t5099) * (-(qJ(3,1) + rSges(3,3)) * t5281 + t5137);
t4695 = 0.1e1 / (((-t4981 * t5156 - t4934) * t5029 + t5030 * t4918) * t5034 + t5144);
t5266 = t4605 * t4695;
t5265 = t4606 * t4702;
t5264 = t4607 * t4703;
t5263 = t4608 * t4704;
t5262 = t4609 * t4705;
t5261 = t4610 * t4706;
t5260 = t4611 * t4707;
t4912 = t4938 * t5051;
t4855 = t5048 * t5193 + t4912;
t4928 = -t4983 * t5027 + t4995 * t5032;
t5118 = t4954 * t5029 * t5030;
t5096 = t4852 * t5035 + t5118;
t5155 = t5035 * t5051;
t5199 = t4954 * t5048;
t5259 = (((t4928 * t5199 + (t4954 * t5155 + t4976 * t5030) * t4922) * t5013 - t5300 * t5007) * t5025 + ((t4855 * t4928 - t4922 * t5096) * t5013 + t5143 * t5007) * t5034 - t4976 * (t5035 * t5007 + t5013 * t5321)) * t5338;
t4913 = t4938 * t5052;
t4856 = t5049 * t5192 + t4913;
t4929 = -t4984 * t5027 + t4996 * t5032;
t5095 = t4853 * t5035 + t5118;
t5154 = t5035 * t5052;
t5198 = t4954 * t5049;
t5258 = (((t4929 * t5198 + (t4954 * t5154 + t4977 * t5030) * t4923) * t5014 - t5301 * t5008) * t5025 + ((t4856 * t4929 - t4923 * t5095) * t5014 + t5142 * t5008) * t5034 - t4977 * (t5008 * t5035 + t5014 * t5320)) * t5339;
t4914 = t4938 * t5053;
t4857 = t5050 * t5191 + t4914;
t4930 = -t4985 * t5027 + t4997 * t5032;
t5094 = t4854 * t5035 + t5118;
t5153 = t5035 * t5053;
t5197 = t4954 * t5050;
t5257 = (((t4930 * t5197 + (t4954 * t5153 + t4978 * t5030) * t4924) * t5015 - t5302 * t5009) * t5025 + ((t4857 * t4930 - t4924 * t5094) * t5015 + t5141 * t5009) * t5034 - t4978 * (t5009 * t5035 + t5015 * t5319)) * t5340;
t4919 = t4938 * t5063;
t4861 = t5060 * t5190 + t4919;
t4931 = -t4986 * t5027 + t4998 * t5032;
t5093 = t4858 * t5035 + t5118;
t5152 = t5035 * t5063;
t5196 = t4954 * t5060;
t5256 = (((t4931 * t5196 + (t4954 * t5152 + t4979 * t5030) * t4925) * t5016 - t5303 * t5010) * t5025 + ((t4861 * t4931 - t4925 * t5093) * t5016 + t5140 * t5010) * t5034 - t4979 * (t5010 * t5035 + t5016 * t5318)) * t5341;
t4920 = t4938 * t5064;
t4862 = t5061 * t5189 + t4920;
t4932 = -t4987 * t5027 + t4999 * t5032;
t5092 = t4859 * t5035 + t5118;
t5151 = t5035 * t5064;
t5195 = t4954 * t5061;
t5255 = (((t4932 * t5195 + (t4954 * t5151 + t4980 * t5030) * t4926) * t5017 - t5304 * t5011) * t5025 + ((t4862 * t4932 - t4926 * t5092) * t5017 + t5139 * t5011) * t5034 - t4980 * (t5011 * t5035 + t5017 * t5317)) * t5342;
t4921 = t4938 * t5065;
t4863 = t5062 * t5188 + t4921;
t4933 = -t4988 * t5027 + t5000 * t5032;
t5091 = t4860 * t5035 + t5118;
t5150 = t5035 * t5065;
t5194 = t4954 * t5062;
t5254 = (((t4933 * t5194 + (t4954 * t5150 + t4981 * t5030) * t4927) * t5018 - t5305 * t5012) * t5025 + ((t4863 * t4933 - t4927 * t5091) * t5018 + t5138 * t5012) * t5034 - t4981 * (t5035 * t5012 + t5018 * t5316)) * t5343;
t4870 = t5020 * t5073 - t5023 * t5299;
t5253 = t5338 * t4870;
t5252 = t5338 * t4882;
t5251 = t4680 * t4870;
t4871 = t5020 * t5074 - t5023 * t5298;
t5250 = t5339 * t4871;
t5249 = t5339 * t4883;
t5248 = t4683 * t4871;
t4872 = t5020 * t5075 - t5023 * t5297;
t5247 = t5340 * t4872;
t5246 = t5340 * t4884;
t5245 = t4686 * t4872;
t4873 = t5020 * t5076 - t5023 * t5296;
t5244 = t5341 * t4873;
t5243 = t5341 * t4885;
t5242 = t4689 * t4873;
t4874 = t5020 * t5077 - t5023 * t5295;
t5241 = t5342 * t4874;
t5240 = t5342 * t4886;
t5239 = t4692 * t4874;
t4875 = t5020 * t5078 - t5023 * t5294;
t5238 = t5343 * t4875;
t5237 = t5343 * t4887;
t5236 = t4695 * t4875;
t5235 = t4702 * t4870;
t5234 = t4703 * t4871;
t5233 = t4704 * t4872;
t5232 = t4705 * t4873;
t5231 = t4706 * t4874;
t5230 = t4707 * t4875;
t5229 = (-t4876 * t5048 + t4912) * t5035;
t5228 = (-t4877 * t5049 + t4913) * t5035;
t5227 = (-t4878 * t5050 + t4914) * t5035;
t5226 = (-t4879 * t5060 + t4919) * t5035;
t5225 = (-t4880 * t5061 + t4920) * t5035;
t5224 = (-t4881 * t5062 + t4921) * t5035;
t5026 = m(1) + m(2) + m(3);
t5223 = t4882 * t5026;
t5222 = t4883 * t5026;
t5221 = t4884 * t5026;
t5220 = t4885 * t5026;
t5219 = t4886 * t5026;
t5218 = t4887 * t5026;
t5217 = t4928 * t5030;
t5216 = t4929 * t5030;
t5215 = t4930 * t5030;
t5214 = t4931 * t5030;
t5213 = t4932 * t5030;
t5212 = t4933 * t5030;
t5174 = t5020 * t5071;
t5173 = t5020 * t5072;
t5166 = t5023 * t5070;
t5165 = t5027 * t5030;
t5164 = t5027 * t5035;
t5163 = t5030 * t5032;
t5130 = (t4928 * (t5034 * t5199 + t4855) * t5035 - t5349 * t4922) * t4702 * t5013;
t5129 = (t4929 * (t5034 * t5198 + t4856) * t5035 - t5348 * t4923) * t4703 * t5014;
t5128 = (t4930 * (t5034 * t5197 + t4857) * t5035 - t5347 * t4924) * t4704 * t5015;
t5127 = (t4931 * (t5034 * t5196 + t4861) * t5035 - t5346 * t4925) * t4705 * t5016;
t5126 = (t4932 * (t5034 * t5195 + t4862) * t5035 - t5345 * t4926) * t4706 * t5017;
t5125 = (t4933 * (t5034 * t5194 + t4863) * t5035 - t5344 * t4927) * t4707 * t5018;
t5124 = t5338 * (t4928 * t5162 + t5029 * (-t4922 * t5048 + t4928 * t5155)) * t5013;
t5123 = t5339 * (t4929 * t5162 + t5029 * (-t4923 * t5049 + t4929 * t5154)) * t5014;
t5122 = t5340 * (t4930 * t5162 + t5029 * (-t4924 * t5050 + t4930 * t5153)) * t5015;
t5121 = t5341 * (t4931 * t5162 + t5029 * (-t4925 * t5060 + t4931 * t5152)) * t5016;
t5120 = t5342 * (t4932 * t5162 + t5029 * (-t4926 * t5061 + t4932 * t5151)) * t5017;
t5119 = t5343 * (t4933 * t5162 + t5029 * (-t4927 * t5062 + t4933 * t5150)) * t5018;
t4906 = t5027 * t5118;
t4839 = -t5012 * t5316 + t5018 * t5035;
t4838 = -t5011 * t5317 + t5017 * t5035;
t4837 = -t5010 * t5318 + t5016 * t5035;
t4836 = -t5009 * t5319 + t5015 * t5035;
t4835 = -t5008 * t5320 + t5014 * t5035;
t4834 = -t5007 * t5321 + t5013 * t5035;
t4827 = t4927 * t5182 - t4933 * t5006;
t4826 = t4926 * t5183 - t4932 * t5005;
t4825 = t4925 * t5184 - t4931 * t5004;
t4824 = t4924 * t5185 - t4930 * t5003;
t4823 = t4923 * t5186 - t4929 * t5002;
t4822 = t4922 * t5187 - t4928 * t5001;
t4821 = t4927 * t5006 + t4933 * t5182;
t4820 = t4926 * t5005 + t4932 * t5183;
t4819 = t4925 * t5004 + t4931 * t5184;
t4818 = t4924 * t5003 + t4930 * t5185;
t4817 = t4923 * t5002 + t4929 * t5186;
t4816 = t4922 * t5001 + t4928 * t5187;
t4815 = -t4927 * t4994 + t4933 * t5176;
t4814 = t4927 * t5176 + t4933 * t4994;
t4813 = -t4926 * t4993 + t4932 * t5177;
t4812 = t4926 * t5177 + t4932 * t4993;
t4811 = -t4925 * t4992 + t4931 * t5178;
t4810 = t4925 * t5178 + t4931 * t4992;
t4809 = -t4924 * t4991 + t4930 * t5179;
t4808 = t4924 * t5179 + t4930 * t4991;
t4807 = -t4923 * t4990 + t4929 * t5180;
t4806 = t4923 * t5180 + t4929 * t4990;
t4805 = -t4922 * t4989 + t4928 * t5181;
t4804 = t4922 * t5181 + t4928 * t4989;
t4755 = -t4981 * t5163 + (t5027 * t5062 - t5032 * t5150) * t4954;
t4754 = -t4981 * t5165 + (-t5027 * t5150 - t5032 * t5062) * t4954;
t4753 = -t4980 * t5163 + (t5027 * t5061 - t5032 * t5151) * t4954;
t4752 = -t4980 * t5165 + (-t5027 * t5151 - t5032 * t5061) * t4954;
t4751 = -t4979 * t5163 + (t5027 * t5060 - t5032 * t5152) * t4954;
t4750 = -t4979 * t5165 + (-t5027 * t5152 - t5032 * t5060) * t4954;
t4749 = -t4978 * t5163 + (t5027 * t5050 - t5032 * t5153) * t4954;
t4748 = -t4978 * t5165 + (-t5027 * t5153 - t5032 * t5050) * t4954;
t4747 = -t4977 * t5163 + (t5027 * t5049 - t5032 * t5154) * t4954;
t4746 = -t4977 * t5165 + (-t5027 * t5154 - t5032 * t5049) * t4954;
t4745 = -t4976 * t5163 + (t5027 * t5048 - t5032 * t5155) * t4954;
t4744 = -t4976 * t5165 + (-t5027 * t5155 - t5032 * t5048) * t4954;
t4719 = t4863 * t5027 + t5032 * t5091;
t4718 = t5027 * t4862 + t5032 * t5092;
t4717 = t5027 * t4861 + t5032 * t5093;
t4716 = -t4860 * t5164 + t4863 * t5032 - t4906;
t4715 = -t4859 * t5164 + t4862 * t5032 - t4906;
t4714 = -t4858 * t5164 + t4861 * t5032 - t4906;
t4713 = t5027 * t4857 + t5032 * t5094;
t4712 = t5027 * t4856 + t5032 * t5095;
t4711 = t5027 * t4855 + t5032 * t5096;
t4710 = -t4854 * t5164 + t4857 * t5032 - t4906;
t4709 = -t4853 * t5164 + t4856 * t5032 - t4906;
t4708 = -t4852 * t5164 + t4855 * t5032 - t4906;
t4677 = -t4754 * t4988 + t4755 * t5000;
t4676 = -t4752 * t4987 + t4753 * t4999;
t4675 = -t4750 * t4986 + t4751 * t4998;
t4674 = -t4748 * t4985 + t4749 * t4997;
t4673 = -t4746 * t4984 + t4747 * t4996;
t4672 = -t4744 * t4983 + t4745 * t4995;
t4659 = -t5305 * t5018 + (t4754 * t5000 + t4755 * t4988) * t5012;
t4658 = -t5304 * t5017 + (t4752 * t4999 + t4753 * t4987) * t5011;
t4657 = -t5303 * t5016 + (t4750 * t4998 + t4751 * t4986) * t5010;
t4656 = -t5302 * t5015 + (t4748 * t4997 + t4749 * t4985) * t5009;
t4655 = -t5301 * t5014 + (t4746 * t4996 + t4747 * t4984) * t5008;
t4654 = -t5300 * t5013 + (t4744 * t4995 + t4745 * t4983) * t5007;
t4653 = -t4821 * t5162 - t5029 * (t4821 * t5150 - t4827 * t5062);
t4652 = -t4820 * t5162 - t5029 * (t4820 * t5151 - t4826 * t5061);
t4651 = -t4819 * t5162 - t5029 * (t4819 * t5152 - t4825 * t5060);
t4650 = t4815 * t5162 + t5029 * (-t4814 * t5062 + t4815 * t5150);
t4649 = t4813 * t5162 + t5029 * (-t4812 * t5061 + t4813 * t5151);
t4648 = t4811 * t5162 + t5029 * (-t4810 * t5060 + t4811 * t5152);
t4647 = -t4818 * t5162 - t5029 * (t4818 * t5153 - t4824 * t5050);
t4646 = -t4817 * t5162 - t5029 * (t4817 * t5154 - t4823 * t5049);
t4645 = -t4816 * t5162 - t5029 * (t4816 * t5155 - t4822 * t5048);
t4644 = t4809 * t5162 + t5029 * (-t4808 * t5050 + t4809 * t5153);
t4643 = t4807 * t5162 + t5029 * (-t4806 * t5049 + t4807 * t5154);
t4642 = t4805 * t5162 + t5029 * (-t4804 * t5048 + t4805 * t5155);
t4641 = -t4821 * t5224 + t4827 * t5344;
t4640 = -t4814 * t5344 + t4815 * t5224;
t4639 = -t4820 * t5225 + t4826 * t5345;
t4638 = -t4812 * t5345 + t4813 * t5225;
t4637 = -t4819 * t5226 + t4825 * t5346;
t4636 = -t4810 * t5346 + t4811 * t5226;
t4635 = t4716 * t4988 + t4719 * t5000;
t4634 = t4715 * t4987 + t4718 * t4999;
t4633 = t4714 * t4986 + t4717 * t4998;
t4632 = -t4818 * t5227 + t4824 * t5347;
t4631 = -t4808 * t5347 + t4809 * t5227;
t4630 = -t4817 * t5228 + t4823 * t5348;
t4629 = -t4806 * t5348 + t4807 * t5228;
t4628 = -t4816 * t5229 + t4822 * t5349;
t4627 = -t4804 * t5349 + t4805 * t5229;
t4626 = t4710 * t4985 + t4713 * t4997;
t4625 = t4709 * t4984 + t4712 * t4996;
t4624 = t4708 * t4983 + t4711 * t4995;
t4623 = -t5138 * t5018 + (t4716 * t5000 - t4719 * t4988) * t5012;
t4622 = -t5139 * t5017 + (t4715 * t4999 - t4718 * t4987) * t5011;
t4621 = -t5140 * t5016 + (t4714 * t4998 - t4717 * t4986) * t5010;
t4620 = -t5141 * t5015 + (t4710 * t4997 - t4713 * t4985) * t5009;
t4619 = -t5142 * t5014 + (t4709 * t4996 - t4712 * t4984) * t5008;
t4618 = -t5143 * t5013 + (t4708 * t4995 - t4711 * t4983) * t5007;
t4599 = (-t4659 * t4994 + t4677 * t5006) * t5025 + (t4623 * t4994 + t4635 * t5006) * t5034 + t4981 * (t4839 * t4994 + t5006 * t5212);
t4598 = (t4659 * t5006 + t4677 * t4994) * t5025 + (-t4623 * t5006 + t4635 * t4994) * t5034 - t4981 * (t4839 * t5006 - t4994 * t5212);
t4597 = (-t4658 * t4993 + t4676 * t5005) * t5025 + (t4622 * t4993 + t4634 * t5005) * t5034 + t4980 * (t4838 * t4993 + t5005 * t5213);
t4596 = (t4658 * t5005 + t4676 * t4993) * t5025 + (-t4622 * t5005 + t4634 * t4993) * t5034 - t4980 * (t4838 * t5005 - t4993 * t5213);
t4595 = (-t4657 * t4992 + t4675 * t5004) * t5025 + (t4621 * t4992 + t4633 * t5004) * t5034 + t4979 * (t4837 * t4992 + t5004 * t5214);
t4594 = (t4657 * t5004 + t4675 * t4992) * t5025 + (-t4621 * t5004 + t4633 * t4992) * t5034 - t4979 * (t4837 * t5004 - t4992 * t5214);
t4593 = (-t4656 * t4991 + t4674 * t5003) * t5025 + (t4620 * t4991 + t4626 * t5003) * t5034 + t4978 * (t4836 * t4991 + t5003 * t5215);
t4592 = (t4656 * t5003 + t4674 * t4991) * t5025 + (-t4620 * t5003 + t4626 * t4991) * t5034 - t4978 * (t4836 * t5003 - t4991 * t5215);
t4591 = (-t4655 * t4990 + t4673 * t5002) * t5025 + (t4619 * t4990 + t4625 * t5002) * t5034 + t4977 * (t4835 * t4990 + t5002 * t5216);
t4590 = (t4655 * t5002 + t4673 * t4990) * t5025 + (-t4619 * t5002 + t4625 * t4990) * t5034 - t4977 * (t4835 * t5002 - t4990 * t5216);
t4589 = (-t4654 * t4989 + t4672 * t5001) * t5025 + (t4618 * t4989 + t4624 * t5001) * t5034 + t4976 * (t4834 * t4989 + t5001 * t5217);
t4588 = (t4654 * t5001 + t4672 * t4989) * t5025 + (-t4618 * t5001 + t4624 * t4989) * t5034 - t4976 * (t4834 * t5001 - t4989 * t5217);
t1 = [-t4600 * t5124 - t4601 * t5123 - t4602 * t5122 - t4603 * t5121 - t4604 * t5120 - t4605 * t5119 - m(4) * g(1) + (-t4882 * t5259 - t4883 * t5258 - t4884 * t5257 - t4885 * t5256 - t4886 * t5255 - t4887 * t5254) * t5026 + (t4606 * t5130 + t4607 * t5129 + t4608 * t5128 + t4609 * t5127 + t4610 * t5126 + t4611 * t5125) * m(3); t4645 * t5271 + t4646 * t5270 + t4647 * t5269 + t4651 * t5268 + t4652 * t5267 + t4653 * t5266 - m(4) * g(2) + (-t4589 * t5252 - t4591 * t5249 - t4593 * t5246 - t4595 * t5243 - t4597 * t5240 - t4599 * t5237) * t5026 + (-t4628 * t5265 - t4630 * t5264 - t4632 * t5263 - t4637 * t5262 - t4639 * t5261 - t4641 * t5260) * m(3); t4642 * t5271 + t4643 * t5270 + t4644 * t5269 + t4648 * t5268 + t4649 * t5267 + t4650 * t5266 - m(4) * g(3) + (-t4588 * t5252 - t4590 * t5249 - t4592 * t5246 - t4594 * t5243 - t4596 * t5240 - t4598 * t5237) * t5026 + (-t4627 * t5265 - t4629 * t5264 - t4631 * t5263 - t4636 * t5262 - t4638 * t5261 - t4640 * t5260) * m(3); -(-t4598 * t4800 - t4599 * t4801) * t5343 * t5218 + (-t4650 * t4800 - t4653 * t4801) * t5266 - (-t4596 * t4796 - t4597 * t4797) * t5342 * t5219 + (-t4649 * t4796 - t4652 * t4797) * t5267 - (-t4594 * t4792 - t4595 * t4793) * t5341 * t5220 + (-t4648 * t4792 - t4651 * t4793) * t5268 - (-t4592 * t4788 - t4593 * t4789) * t5340 * t5221 + (-t4644 * t4788 - t4647 * t4789) * t5269 - (-t4590 * t4784 - t4591 * t4785) * t5339 * t5222 + (-t4643 * t4784 - t4646 * t4785) * t5270 - (-t4588 * t4780 - t4589 * t4781) * t5338 * t5223 + (-t4642 * t4780 - t4645 * t4781) * t5271 + (((-g(2) * t5173 - g(3) * t5071) * t5022 + g(2) * t5166) * t5024 + ((g(2) * t5071 - g(3) * t5173) * t5022 + g(3) * t5166) * t5021 + ((g(2) * t5174 - g(3) * t5072) * t5024 + (g(2) * t5072 + g(3) * t5174) * t5021) * t5019) * m(4) + (-(-t4640 * t4800 - t4641 * t4801) * t5260 - (-t4638 * t4796 - t4639 * t4797) * t5261 - (-t4636 * t4792 - t4637 * t4793) * t5262 - (-t4631 * t4788 - t4632 * t4789) * t5263 - (-t4629 * t4784 - t4630 * t4785) * t5264 - (-t4627 * t4780 - t4628 * t4781) * t5265) * m(3); -(-t4598 * t5238 + t4801 * t5254) * t5218 + (-t4650 * t5236 - t4801 * t5119) * t4605 - (-t4640 * t5230 - t4801 * t5125) * t5282 - (-t4596 * t5241 + t4797 * t5255) * t5219 + (-t4649 * t5239 - t4797 * t5120) * t4604 - (-t4638 * t5231 - t4797 * t5126) * t5283 - (-t4594 * t5244 + t4793 * t5256) * t5220 + (-t4648 * t5242 - t4793 * t5121) * t4603 - (-t4636 * t5232 - t4793 * t5127) * t5284 - (-t4592 * t5247 + t4789 * t5257) * t5221 + (-t4644 * t5245 - t4789 * t5122) * t4602 - (-t4631 * t5233 - t4789 * t5128) * t5285 - (-t4590 * t5250 + t4785 * t5258) * t5222 + (-t4643 * t5248 - t4785 * t5123) * t4601 - (-t4629 * t5234 - t4785 * t5129) * t5286 - (-t4588 * t5253 + t4781 * t5259) * t5223 + (-t4642 * t5251 - t4781 * t5124) * t4600 - (-t4627 * t5235 - t4781 * t5130) * t5287 - (t5331 * g(3) + ((t5020 * t5098 + t5166) * t5024 + (t5019 * t5072 + t5022 * t5071) * t5021) * g(1)) * m(4); -(t4599 * t5238 + t4800 * t5254) * t5218 + (t4653 * t5236 - t4800 * t5119) * t4605 - (t4641 * t5230 - t4800 * t5125) * t5282 - (t4597 * t5241 + t4796 * t5255) * t5219 + (t4652 * t5239 - t4796 * t5120) * t4604 - (t4639 * t5231 - t4796 * t5126) * t5283 - (t4595 * t5244 + t4792 * t5256) * t5220 + (t4651 * t5242 - t4792 * t5121) * t4603 - (t4637 * t5232 - t4792 * t5127) * t5284 - (t4593 * t5247 + t4788 * t5257) * t5221 + (t4647 * t5245 - t4788 * t5122) * t4602 - (t4632 * t5233 - t4788 * t5128) * t5285 - (t4591 * t5250 + t4784 * t5258) * t5222 + (t4646 * t5248 - t4784 * t5123) * t4601 - (t4630 * t5234 - t4784 * t5129) * t5286 - (t4589 * t5253 + t4780 * t5259) * t5223 + (t4645 * t5251 - t4780 * t5124) * t4600 - (t4628 * t5235 - t4780 * t5130) * t5287 + (t5331 * g(2) + (-t5021 * t5166 + (t5021 * t5173 + t5024 * t5071) * t5022 + (-t5021 * t5174 + t5024 * t5072) * t5019) * g(1)) * m(4);];
taugX  = t1;
