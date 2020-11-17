% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:13:05
% EndTime: 2020-08-07 11:13:11
% DurationCPUTime: 5.82s
% Computational Cost: add. (2210->425), mult. (4666->810), div. (112->17), fcn. (4976->30), ass. (0->301)
t5111 = MDP(1) * g(3);
t4981 = cos(qJ(2,1));
t4975 = sin(qJ(2,1));
t4982 = pkin(7) + pkin(6);
t5028 = t4975 * t4982;
t4929 = pkin(2) * t4981 + t5028;
t4958 = sin(pkin(8));
t4960 = cos(pkin(8));
t4939 = t4982 * t4981;
t4926 = pkin(2) * t4975 - t4939;
t4961 = cos(pkin(4));
t4959 = sin(pkin(4));
t4974 = sin(qJ(3,1));
t5062 = t4959 * t4974;
t4993 = pkin(3) * t5062 - t4926 * t4961;
t5115 = t4929 * t4960 + t4993 * t4958;
t4979 = cos(qJ(2,2));
t4973 = sin(qJ(2,2));
t5032 = t4973 * t4982;
t4928 = pkin(2) * t4979 + t5032;
t4938 = t4982 * t4979;
t4925 = pkin(2) * t4973 - t4938;
t4972 = sin(qJ(3,2));
t5064 = t4959 * t4972;
t4994 = pkin(3) * t5064 - t4925 * t4961;
t5114 = t4928 * t4960 + t4994 * t4958;
t4977 = cos(qJ(2,3));
t4971 = sin(qJ(2,3));
t5036 = t4971 * t4982;
t4927 = pkin(2) * t4977 + t5036;
t4937 = t4982 * t4977;
t4924 = pkin(2) * t4971 - t4937;
t4970 = sin(qJ(3,3));
t5066 = t4959 * t4970;
t4995 = pkin(3) * t5066 - t4924 * t4961;
t5113 = t4927 * t4960 + t4995 * t4958;
t4969 = cos(qJ(2,4));
t4967 = sin(qJ(2,4));
t5041 = t4967 * t4982;
t4923 = pkin(2) * t4969 + t5041;
t4933 = t4982 * t4969;
t4922 = pkin(2) * t4967 - t4933;
t4966 = sin(qJ(3,4));
t5069 = t4959 * t4966;
t4996 = pkin(3) * t5069 - t4922 * t4961;
t5112 = t4923 * t4960 + t4996 * t4958;
t4968 = cos(qJ(3,4));
t5110 = pkin(3) * t4968 ^ 2;
t4976 = cos(qJ(3,3));
t5109 = pkin(3) * t4976 ^ 2;
t4978 = cos(qJ(3,2));
t5108 = pkin(3) * t4978 ^ 2;
t4980 = cos(qJ(3,1));
t5107 = pkin(3) * t4980 ^ 2;
t5106 = g(3) * t4959;
t4984 = koppelP(4,2);
t4988 = koppelP(4,1);
t4912 = t4958 * t4984 + t4960 * t4988;
t4913 = -t4958 * t4988 + t4960 * t4984;
t4983 = xP(4);
t4952 = sin(t4983);
t4953 = cos(t4983);
t4852 = t4912 * t4952 + t4913 * t4953;
t4962 = legFrame(4,3);
t4940 = sin(t4962);
t4944 = cos(t4962);
t5000 = t4912 * t4953 - t4913 * t4952;
t4796 = t4852 * t4944 - t5000 * t4940;
t5056 = t4961 * t4967;
t4892 = t4969 * t4984 - t4988 * t5056;
t4893 = t4969 * t4988 + t4984 * t5056;
t4836 = t4892 * t4960 - t4893 * t4958;
t4837 = t4892 * t4958 + t4893 * t4960;
t4932 = pkin(3) * t4968 + pkin(2);
t5057 = t4961 * t4966;
t4908 = t4932 * t5057;
t5042 = t4967 * t4968;
t5067 = t4959 * t4968;
t5105 = (((-t4836 * t4952 + t4837 * t4953) * t4944 + (t4836 * t4953 + t4837 * t4952) * t4940) * t4966 + t4796 * t5067) / (t4908 + (pkin(3) * t5042 + t4922) * t5067);
t4985 = koppelP(3,2);
t4989 = koppelP(3,1);
t4914 = t4958 * t4985 + t4960 * t4989;
t4915 = -t4958 * t4989 + t4960 * t4985;
t4853 = t4914 * t4952 + t4915 * t4953;
t4963 = legFrame(3,3);
t4941 = sin(t4963);
t4945 = cos(t4963);
t4999 = t4914 * t4953 - t4915 * t4952;
t4797 = t4853 * t4945 - t4999 * t4941;
t5053 = t4961 * t4971;
t4894 = t4977 * t4985 - t4989 * t5053;
t4897 = t4977 * t4989 + t4985 * t5053;
t4838 = t4894 * t4960 - t4897 * t4958;
t4841 = t4894 * t4958 + t4897 * t4960;
t4934 = pkin(3) * t4976 + pkin(2);
t5054 = t4961 * t4970;
t4909 = t4934 * t5054;
t5037 = t4971 * t4976;
t5060 = t4959 * t4976;
t5104 = (((-t4838 * t4952 + t4841 * t4953) * t4945 + (t4838 * t4953 + t4841 * t4952) * t4941) * t4970 + t4797 * t5060) / (t4909 + (pkin(3) * t5037 + t4924) * t5060);
t4986 = koppelP(2,2);
t4990 = koppelP(2,1);
t4916 = t4958 * t4986 + t4960 * t4990;
t4917 = -t4958 * t4990 + t4960 * t4986;
t4854 = t4916 * t4952 + t4917 * t4953;
t4964 = legFrame(2,3);
t4942 = sin(t4964);
t4946 = cos(t4964);
t4998 = t4916 * t4953 - t4917 * t4952;
t4798 = t4854 * t4946 - t4998 * t4942;
t5051 = t4961 * t4973;
t4895 = t4979 * t4986 - t4990 * t5051;
t4898 = t4979 * t4990 + t4986 * t5051;
t4839 = t4895 * t4960 - t4898 * t4958;
t4842 = t4895 * t4958 + t4898 * t4960;
t4935 = pkin(3) * t4978 + pkin(2);
t5052 = t4961 * t4972;
t4910 = t4935 * t5052;
t5033 = t4973 * t4978;
t5059 = t4959 * t4978;
t5103 = (((-t4839 * t4952 + t4842 * t4953) * t4946 + t4942 * (t4839 * t4953 + t4842 * t4952)) * t4972 + t4798 * t5059) / (t4910 + (pkin(3) * t5033 + t4925) * t5059);
t4987 = koppelP(1,2);
t4991 = koppelP(1,1);
t4918 = t4958 * t4987 + t4960 * t4991;
t4919 = -t4958 * t4991 + t4960 * t4987;
t4855 = t4918 * t4952 + t4919 * t4953;
t4965 = legFrame(1,3);
t4943 = sin(t4965);
t4947 = cos(t4965);
t4997 = t4918 * t4953 - t4919 * t4952;
t4799 = t4855 * t4947 - t4943 * t4997;
t5049 = t4961 * t4975;
t4896 = t4981 * t4987 - t4991 * t5049;
t4899 = t4981 * t4991 + t4987 * t5049;
t4840 = t4896 * t4960 - t4899 * t4958;
t4843 = t4896 * t4958 + t4899 * t4960;
t4936 = pkin(3) * t4980 + pkin(2);
t5050 = t4961 * t4974;
t4911 = t4936 * t5050;
t5029 = t4975 * t4980;
t5058 = t4959 * t4980;
t5102 = (((-t4840 * t4952 + t4843 * t4953) * t4947 + (t4840 * t4953 + t4843 * t4952) * t4943) * t4974 + t4799 * t5058) / (t4911 + (pkin(3) * t5029 + t4926) * t5058);
t4876 = t4932 * t4967 - t4933;
t4844 = 0.1e1 / (t4876 * t5067 + t4908);
t5077 = (t4932 * t4969 + t5041) * t4961;
t5101 = (t4796 * t5077 - t4876 * (t4852 * t4940 + t5000 * t4944)) * t4844;
t4889 = t4934 * t4971 - t4937;
t4846 = 0.1e1 / (t4889 * t5060 + t4909);
t5076 = (t4934 * t4977 + t5036) * t4961;
t5100 = (t4797 * t5076 - t4889 * (t4853 * t4941 + t4999 * t4945)) * t4846;
t4890 = t4935 * t4973 - t4938;
t4847 = 0.1e1 / (t4890 * t5059 + t4910);
t5075 = (t4935 * t4979 + t5032) * t4961;
t5099 = (t4798 * t5075 - (t4854 * t4942 + t4998 * t4946) * t4890) * t4847;
t4891 = t4936 * t4975 - t4939;
t4848 = 0.1e1 / (t4891 * t5058 + t4911);
t5074 = (t4936 * t4981 + t5028) * t4961;
t5098 = (t4799 * t5074 - (t4943 * t4855 + t4997 * t4947) * t4891) * t4848;
t4930 = g(1) * t4958 - g(2) * t4960;
t4931 = g(1) * t4960 + g(2) * t4958;
t4792 = (-t5106 + (t4930 * t4944 + t4931 * t4940) * t4961) * t4969 + (-t4930 * t4940 + t4931 * t4944) * t4967;
t5068 = t4959 * t4967;
t4824 = 0.1e1 / (t5068 * t5110 + (pkin(3) * t5057 + t4922 * t4959) * t4968 + pkin(2) * t5057);
t5097 = t4792 * t4824;
t4793 = (-t5106 + (t4930 * t4945 + t4931 * t4941) * t4961) * t4977 + (-t4930 * t4941 + t4931 * t4945) * t4971;
t5065 = t4959 * t4971;
t4825 = 0.1e1 / (t5065 * t5109 + (pkin(3) * t5054 + t4924 * t4959) * t4976 + pkin(2) * t5054);
t5096 = t4793 * t4825;
t4794 = (-t5106 + (t4930 * t4946 + t4931 * t4942) * t4961) * t4979 + (-t4930 * t4942 + t4931 * t4946) * t4973;
t5063 = t4959 * t4973;
t4826 = 0.1e1 / (t5063 * t5108 + (pkin(3) * t5052 + t4925 * t4959) * t4978 + pkin(2) * t5052);
t5095 = t4794 * t4826;
t4795 = (-t5106 + (t4930 * t4947 + t4931 * t4943) * t4961) * t4981 + (-t4930 * t4943 + t4931 * t4947) * t4975;
t5061 = t4959 * t4975;
t4827 = 0.1e1 / (t5061 * t5107 + (pkin(3) * t5050 + t4926 * t4959) * t4980 + pkin(2) * t5050);
t5094 = t4795 * t4827;
t4856 = -t4940 * t4958 + t4944 * t4960;
t4860 = t4940 * t4960 + t4944 * t4958;
t5093 = (-t4856 * t5077 + t4860 * t4876) * t4844;
t5092 = (-t4856 * t4876 - t4860 * t5077) * t4844;
t4857 = -t4941 * t4958 + t4945 * t4960;
t4861 = t4941 * t4960 + t4945 * t4958;
t5091 = (-t4857 * t5076 + t4861 * t4889) * t4846;
t4858 = -t4942 * t4958 + t4946 * t4960;
t4862 = t4942 * t4960 + t4946 * t4958;
t5090 = (-t4858 * t5075 + t4862 * t4890) * t4847;
t4859 = -t4943 * t4958 + t4947 * t4960;
t4863 = t4943 * t4960 + t4947 * t4958;
t5089 = (-t4859 * t5074 + t4863 * t4891) * t4848;
t5088 = (-t4857 * t4889 - t4861 * t5076) * t4846;
t5087 = (-t4858 * t4890 - t4862 * t5075) * t4847;
t5086 = (-t4859 * t4891 - t4863 * t5074) * t4848;
t4900 = g(1) * t4940 - g(2) * t4944;
t4904 = g(1) * t4944 + g(2) * t4940;
t5055 = t4961 * t4969;
t4816 = t4904 * (t4958 * t5055 + t4960 * t4967) + t4900 * (-t4958 * t4967 + t4960 * t5055) - t4969 * t5106;
t5085 = t4816 * t4824;
t4864 = t4958 * t5056 - t4960 * t4969;
t4865 = t4958 * t4969 + t4960 * t5056;
t4817 = g(3) * t5068 - t4864 * t4904 - t4865 * t4900;
t5084 = t4817 * t4824;
t4901 = g(1) * t4941 - g(2) * t4945;
t4905 = g(1) * t4945 + g(2) * t4941;
t5048 = t4961 * t4977;
t4818 = t4905 * (t4958 * t5048 + t4960 * t4971) + t4901 * (-t4958 * t4971 + t4960 * t5048) - t4977 * t5106;
t5083 = t4818 * t4825;
t4869 = t4958 * t5053 - t4960 * t4977;
t4872 = t4958 * t4977 + t4960 * t5053;
t4819 = g(3) * t5065 - t4869 * t4905 - t4872 * t4901;
t5082 = t4819 * t4825;
t4902 = g(1) * t4942 - g(2) * t4946;
t4906 = g(1) * t4946 + g(2) * t4942;
t5047 = t4961 * t4979;
t4820 = t4906 * (t4958 * t5047 + t4960 * t4973) + t4902 * (-t4958 * t4973 + t4960 * t5047) - t4979 * t5106;
t5081 = t4820 * t4826;
t4870 = t4958 * t5051 - t4960 * t4979;
t4873 = t4958 * t4979 + t4960 * t5051;
t4821 = g(3) * t5063 - t4870 * t4906 - t4873 * t4902;
t5080 = t4821 * t4826;
t4903 = g(1) * t4943 - g(2) * t4947;
t4907 = g(1) * t4947 + g(2) * t4943;
t5046 = t4961 * t4981;
t4822 = t4907 * (t4958 * t5046 + t4960 * t4975) + t4903 * (-t4958 * t4975 + t4960 * t5046) - t4981 * t5106;
t5079 = t4822 * t4827;
t4871 = t4958 * t5049 - t4960 * t4981;
t4874 = t4958 * t4981 + t4960 * t5049;
t4823 = g(3) * t5061 - t4871 * t4907 - t4874 * t4903;
t5078 = t4823 * t4827;
t5045 = t4961 * t4982;
t5044 = t4966 * t4967;
t5043 = t4966 * t4969;
t5040 = t4968 * t4969;
t5039 = t4970 * t4971;
t5038 = t4970 * t4977;
t5035 = t4972 * t4973;
t5034 = t4972 * t4979;
t5031 = t4974 * t4975;
t5030 = t4974 * t4981;
t5027 = t4976 * t4977;
t5026 = t4978 * t4979;
t5025 = t4980 * t4981;
t5024 = pkin(2) * t5069;
t5023 = pkin(2) * t5066;
t5022 = pkin(2) * t5064;
t5021 = pkin(2) * t5062;
t5020 = t4792 * t5105;
t5019 = t4793 * t5104;
t5018 = t4794 * t5103;
t5017 = t4795 * t5102;
t5016 = t4966 * t5097;
t5015 = t4968 * t5097;
t5014 = t4970 * t5096;
t5013 = t4976 * t5096;
t5012 = t4972 * t5095;
t5011 = t4978 * t5095;
t5010 = t4974 * t5094;
t5009 = t4980 * t5094;
t5008 = t4932 * t5069;
t5007 = t4934 * t5066;
t5006 = t4935 * t5064;
t5005 = t4936 * t5062;
t5004 = t4932 * t4961;
t5003 = t4934 * t4961;
t5002 = t4935 * t4961;
t5001 = t4936 * t4961;
t4992 = 0.1e1 / pkin(3);
t4921 = g(1) * t4953 + g(2) * t4952;
t4920 = g(1) * t4952 - g(2) * t4953;
t4885 = t4961 * t5029 - t5062;
t4884 = t4961 * t5033 - t5064;
t4883 = t4961 * t5037 - t5066;
t4882 = t4959 * t5029 + t5050;
t4881 = t4961 * t5031 + t5058;
t4880 = t4959 * t5033 + t5052;
t4879 = t4961 * t5035 + t5059;
t4878 = t4959 * t5037 + t5054;
t4877 = t4961 * t5039 + t5060;
t4868 = t4961 * t5042 - t5069;
t4867 = t4959 * t5042 + t5057;
t4866 = t4961 * t5044 + t5067;
t4831 = t4929 * t4958 - t4993 * t4960;
t4830 = t4928 * t4958 - t4994 * t4960;
t4829 = t4927 * t4958 - t4995 * t4960;
t4828 = t4923 * t4958 - t4996 * t4960;
t4807 = -t4863 * t5058 - (-t4859 * t4981 + t4863 * t5049) * t4974;
t4806 = -t4862 * t5059 - (-t4858 * t4979 + t4862 * t5051) * t4972;
t4805 = -t4861 * t5060 - (-t4857 * t4977 + t4861 * t5053) * t4970;
t4804 = -t4859 * t5058 - (t4859 * t5049 + t4863 * t4981) * t4974;
t4803 = -t4858 * t5059 - (t4858 * t5051 + t4862 * t4979) * t4972;
t4802 = -t4857 * t5060 - (t4857 * t5053 + t4861 * t4977) * t4970;
t4801 = -t4860 * t5067 - (-t4856 * t4969 + t4860 * t5056) * t4966;
t4800 = -t4856 * t5067 - (t4856 * t5056 + t4860 * t4969) * t4966;
t4791 = (-t4885 * t4958 + t4960 * t5025) * t4907 - t4903 * (t4885 * t4960 + t4958 * t5025) + g(3) * t4882;
t4790 = (-t4884 * t4958 + t4960 * t5026) * t4906 - (t4884 * t4960 + t4958 * t5026) * t4902 + g(3) * t4880;
t4789 = t4907 * (-t4881 * t4958 + t4960 * t5030) - t4903 * (t4881 * t4960 + t4958 * t5030) - g(3) * (-t4959 * t5031 + t4961 * t4980);
t4788 = t4906 * (-t4879 * t4958 + t4960 * t5034) - (t4879 * t4960 + t4958 * t5034) * t4902 - g(3) * (-t4959 * t5035 + t4961 * t4978);
t4787 = t4905 * (-t4883 * t4958 + t4960 * t5027) - (t4883 * t4960 + t4958 * t5027) * t4901 + g(3) * t4878;
t4786 = t4905 * (-t4877 * t4958 + t4960 * t5038) - (t4877 * t4960 + t4958 * t5038) * t4901 - g(3) * (-t4959 * t5039 + t4961 * t4976);
t4785 = t4904 * (-t4868 * t4958 + t4960 * t5040) - (t4868 * t4960 + t4958 * t5040) * t4900 + g(3) * t4867;
t4784 = t4904 * (-t4866 * t4958 + t4960 * t5043) - (t4866 * t4960 + t4958 * t5043) * t4900 - g(3) * (-t4959 * t5044 + t4961 * t4968);
t1 = [(-(-(t4871 * t4947 + t4874 * t4943) * t5107 + (-t4943 * t4831 + t5115 * t4947) * t4980 + t4863 * t5021) * t4827 - (-(t4870 * t4946 + t4873 * t4942) * t5108 + (-t4942 * t4830 + t5114 * t4946) * t4978 + t4862 * t5022) * t4826 - (-(t4869 * t4945 + t4872 * t4941) * t5109 + (-t4941 * t4829 + t5113 * t4945) * t4976 + t4861 * t5023) * t4825 - (-(t4864 * t4944 + t4865 * t4940) * t5110 + (-t4940 * t4828 + t5112 * t4944) * t4968 + t4860 * t5024) * t4824) * t5111 + (t4800 * t5085 + t4802 * t5083 + t4803 * t5081 + t4804 * t5079) * MDP(3) + (t4800 * t5084 + t4802 * t5082 + t4803 * t5080 + t4804 * t5078) * MDP(4) + (t4800 * t5015 + t4802 * t5013 + t4803 * t5011 + t4804 * t5009 + (t4784 * t5093 + t4786 * t5091 + t4788 * t5090 + t4789 * t5089) * t4992) * MDP(10) + (-t4800 * t5016 - t4802 * t5014 - t4803 * t5012 - t4804 * t5010 + (t4785 * t5093 + t4787 * t5091 + t4790 * t5090 + t4791 * t5089) * t4992) * MDP(11) + (-t4920 * t4952 - t4921 * t4953) * MDP(15); (-((-t4871 * t4943 + t4874 * t4947) * t5107 + (t4831 * t4947 + t5115 * t4943) * t4980 - t4859 * t5021) * t4827 - ((-t4870 * t4942 + t4873 * t4946) * t5108 + (t4830 * t4946 + t5114 * t4942) * t4978 - t4858 * t5022) * t4826 - ((-t4869 * t4941 + t4872 * t4945) * t5109 + (t4829 * t4945 + t5113 * t4941) * t4976 - t4857 * t5023) * t4825 - ((-t4864 * t4940 + t4865 * t4944) * t5110 + (t4828 * t4944 + t5112 * t4940) * t4968 - t4856 * t5024) * t4824) * t5111 + (t4801 * t5085 + t4805 * t5083 + t4806 * t5081 + t4807 * t5079) * MDP(3) + (t4801 * t5084 + t4805 * t5082 + t4806 * t5080 + t4807 * t5078) * MDP(4) + (t4801 * t5015 + t4805 * t5013 + t4806 * t5011 + t4807 * t5009 + (t4784 * t5092 + t4786 * t5088 + t4788 * t5087 + t4789 * t5086) * t4992) * MDP(10) + (-t4801 * t5016 - t4805 * t5014 - t4806 * t5012 - t4807 * t5010 + (t4785 * t5092 + t4787 * t5088 + t4790 * t5087 + t4791 * t5086) * t4992) * MDP(11) + (t4920 * t4953 - t4921 * t4952) * MDP(15); (-(4 * MDP(1)) - MDP(15)) * g(3); (-(-(t4952 * t4991 + t4953 * t4987) * ((t4936 * t4859 + t4863 * t5045) * t5025 - (-t4859 * t4982 + t4863 * t5001) * t5029 + t4863 * t5005) + (-t4952 * t4987 + t4953 * t4991) * ((-t4859 * t5045 + t4936 * t4863) * t5025 + (t4859 * t5001 + t4863 * t4982) * t5029 - t4859 * t5005)) / (t4882 * t4936 - t4939 * t5058) - (-(t4952 * t4990 + t4953 * t4986) * ((t4935 * t4858 + t4862 * t5045) * t5026 - (-t4858 * t4982 + t4862 * t5002) * t5033 + t4862 * t5006) + (-t4952 * t4986 + t4953 * t4990) * ((-t4858 * t5045 + t4935 * t4862) * t5026 + (t4858 * t5002 + t4862 * t4982) * t5033 - t4858 * t5006)) / (t4880 * t4935 - t4938 * t5059) - (-(t4952 * t4989 + t4953 * t4985) * ((t4934 * t4857 + t4861 * t5045) * t5027 - (-t4857 * t4982 + t4861 * t5003) * t5037 + t4861 * t5007) + (-t4952 * t4985 + t4953 * t4989) * ((-t4857 * t5045 + t4934 * t4861) * t5027 + (t4857 * t5003 + t4861 * t4982) * t5037 - t4857 * t5007)) / (t4878 * t4934 - t4937 * t5060) - (-(t4952 * t4988 + t4953 * t4984) * ((t4932 * t4856 + t4860 * t5045) * t5040 - (-t4856 * t4982 + t4860 * t5004) * t5042 + t4860 * t5008) + (-t4952 * t4984 + t4953 * t4988) * ((-t4856 * t5045 + t4932 * t4860) * t5040 + (t4856 * t5004 + t4860 * t4982) * t5042 - t4856 * t5008)) / (t4867 * t4932 - t4933 * t5067)) * t5111 + (t4816 * t5105 + t4818 * t5104 + t4820 * t5103 + t4822 * t5102) * MDP(3) + (t4817 * t5105 + t4819 * t5104 + t4821 * t5103 + t4823 * t5102) * MDP(4) + (t4968 * t5020 + t4976 * t5019 + t4978 * t5018 + t4980 * t5017 + (t4784 * t5101 + t4786 * t5100 + t4788 * t5099 + t4789 * t5098) * t4992) * MDP(10) + (-t4966 * t5020 - t4970 * t5019 - t4972 * t5018 - t4974 * t5017 + (t4785 * t5101 + t4787 * t5100 + t4790 * t5099 + t4791 * t5098) * t4992) * MDP(11) + t4920 * MDP(13) + t4921 * MDP(14);];
taugX  = t1;