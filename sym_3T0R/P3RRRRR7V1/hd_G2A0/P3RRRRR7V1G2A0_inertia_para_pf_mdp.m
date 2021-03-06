% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRRRR7V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR7V1G2A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3RRRRR7V1G2A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:50:08
% EndTime: 2020-08-07 03:50:29
% DurationCPUTime: 22.11s
% Computational Cost: add. (18558->799), mult. (32649->1359), div. (3480->27), fcn. (26799->60), ass. (0->580)
t4614 = legFrame(3,2);
t4584 = sin(t4614);
t4587 = cos(t4614);
t4995 = t4584 * t4587;
t4615 = legFrame(2,2);
t4585 = sin(t4615);
t4588 = cos(t4615);
t4994 = t4585 * t4588;
t4616 = legFrame(1,2);
t4586 = sin(t4616);
t4589 = cos(t4616);
t4993 = t4586 * t4589;
t4992 = -2 * pkin(1);
t4991 = 2 * pkin(1);
t4990 = 2 * pkin(2);
t4989 = 2 * MDP(5);
t4988 = 4 * MDP(12);
t4635 = (pkin(4) + pkin(5));
t4987 = 2 * t4635;
t4986 = pkin(4) / 0.2e1;
t4626 = cos(qJ(3,3));
t4985 = pkin(2) * t4626;
t4629 = cos(qJ(3,2));
t4984 = pkin(2) * t4629;
t4632 = cos(qJ(3,1));
t4983 = pkin(2) * t4632;
t4618 = sin(qJ(2,3));
t4982 = pkin(4) * t4618;
t4621 = sin(qJ(2,2));
t4981 = pkin(4) * t4621;
t4624 = sin(qJ(2,1));
t4980 = pkin(4) * t4624;
t4979 = -qJ(3,1) + qJ(1,1);
t4978 = qJ(3,1) + qJ(1,1);
t4977 = -qJ(3,2) + qJ(1,2);
t4976 = qJ(3,2) + qJ(1,2);
t4975 = -qJ(3,3) + qJ(1,3);
t4974 = qJ(3,3) + qJ(1,3);
t4646 = 1 / pkin(1);
t4973 = MDP(6) * t4646;
t4972 = MDP(7) * t4646;
t4971 = MDP(8) / pkin(1) ^ 2;
t4970 = 0.2e1 * qJ(2,1) + qJ(1,1);
t4969 = -0.2e1 * qJ(2,1) + qJ(1,1);
t4968 = 0.2e1 * qJ(2,2) + qJ(1,2);
t4967 = -0.2e1 * qJ(2,2) + qJ(1,2);
t4966 = 0.2e1 * qJ(2,3) + qJ(1,3);
t4965 = -0.2e1 * qJ(2,3) + qJ(1,3);
t4964 = MDP(15) * t4646;
t4628 = cos(qJ(1,3));
t4636 = 0.2e1 * qJ(3,3);
t4495 = (sin(qJ(2,3) - t4975) - sin(qJ(2,3) + t4974)) * t4987 + (-cos(0.2e1 * qJ(3,3) - t4965) - cos(t4636 + t4966) - 0.2e1 * t4628) * pkin(2) + (-cos(qJ(3,3) - t4965) - cos(qJ(3,3) + t4966) - cos(t4975) - cos(t4974)) * pkin(1);
t4611 = qJ(2,3) + qJ(3,3);
t4524 = (-sin(t4636 + qJ(2,3)) + t4618) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t4611)) * pkin(1);
t4963 = t4495 / t4524;
t4631 = cos(qJ(1,2));
t4639 = 0.2e1 * qJ(3,2);
t4496 = (sin(qJ(2,2) - t4977) - sin(qJ(2,2) + t4976)) * t4987 + (-cos(0.2e1 * qJ(3,2) - t4967) - cos(t4639 + t4968) - 0.2e1 * t4631) * pkin(2) + (-cos(qJ(3,2) - t4967) - cos(qJ(3,2) + t4968) - cos(t4977) - cos(t4976)) * pkin(1);
t4612 = qJ(2,2) + qJ(3,2);
t4525 = (-sin(t4639 + qJ(2,2)) + t4621) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t4612)) * pkin(1);
t4962 = t4496 / t4525;
t4634 = cos(qJ(1,1));
t4642 = 0.2e1 * qJ(3,1);
t4497 = (sin(qJ(2,1) - t4979) - sin(qJ(2,1) + t4978)) * t4987 + (-cos(0.2e1 * qJ(3,1) - t4969) - cos(t4642 + t4970) - 0.2e1 * t4634) * pkin(2) + (-cos(qJ(3,1) - t4969) - cos(qJ(3,1) + t4970) - cos(t4979) - cos(t4978)) * pkin(1);
t4613 = qJ(2,1) + qJ(3,1);
t4526 = (-sin(t4642 + qJ(2,1)) + t4624) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t4613)) * pkin(1);
t4961 = t4497 / t4526;
t4619 = sin(qJ(1,3));
t4573 = pkin(1) + t4985;
t4627 = cos(qJ(2,3));
t4617 = sin(qJ(3,3));
t4915 = t4617 * t4618;
t4698 = pkin(2) * t4915 - t4573 * t4627;
t4901 = t4628 * t4635;
t4506 = -t4619 * t4698 - t4901;
t4914 = t4617 * t4627;
t4542 = pkin(2) * t4914 + t4573 * t4618;
t4500 = t4506 * t4584 - t4542 * t4587;
t4590 = 0.1e1 / t4617;
t4960 = t4500 * t4590;
t4501 = -t4506 * t4587 - t4542 * t4584;
t4959 = t4501 * t4590;
t4622 = sin(qJ(1,2));
t4575 = pkin(1) + t4984;
t4630 = cos(qJ(2,2));
t4620 = sin(qJ(3,2));
t4911 = t4620 * t4621;
t4697 = pkin(2) * t4911 - t4575 * t4630;
t4898 = t4631 * t4635;
t4507 = -t4622 * t4697 - t4898;
t4910 = t4620 * t4630;
t4543 = pkin(2) * t4910 + t4575 * t4621;
t4502 = t4507 * t4585 - t4543 * t4588;
t4594 = 0.1e1 / t4620;
t4958 = t4502 * t4594;
t4503 = -t4507 * t4588 - t4543 * t4585;
t4957 = t4503 * t4594;
t4625 = sin(qJ(1,1));
t4577 = pkin(1) + t4983;
t4633 = cos(qJ(2,1));
t4623 = sin(qJ(3,1));
t4907 = t4623 * t4624;
t4696 = pkin(2) * t4907 - t4577 * t4633;
t4895 = t4634 * t4635;
t4508 = -t4625 * t4696 - t4895;
t4906 = t4623 * t4633;
t4544 = pkin(2) * t4906 + t4577 * t4624;
t4504 = t4508 * t4586 - t4544 * t4589;
t4598 = 0.1e1 / t4623;
t4956 = t4504 * t4598;
t4505 = -t4508 * t4589 - t4544 * t4586;
t4955 = t4505 * t4598;
t4512 = -t4619 * t4635 + t4628 * t4698;
t4954 = t4512 * t4590;
t4513 = -t4622 * t4635 + t4631 * t4697;
t4953 = t4513 * t4594;
t4514 = -t4625 * t4635 + t4634 * t4696;
t4952 = t4514 * t4598;
t4536 = 0.1e1 / t4698;
t4951 = t4536 * t4590;
t4950 = 0.1e1 / t4698 ^ 2 / t4617 ^ 2;
t4538 = 0.1e1 / t4697;
t4949 = t4538 * t4594;
t4948 = 0.1e1 / t4697 ^ 2 / t4620 ^ 2;
t4540 = 0.1e1 / t4696;
t4947 = t4540 * t4598;
t4946 = 0.1e1 / t4696 ^ 2 / t4623 ^ 2;
t4566 = pkin(2) * cos(t4611) + t4627 * pkin(1);
t4560 = 0.1e1 / t4566;
t4945 = t4560 * t4619;
t4944 = t4560 * t4628;
t4561 = 0.1e1 / t4566 ^ 2;
t4593 = t4619 ^ 2;
t4943 = t4561 * t4593;
t4604 = t4628 ^ 2;
t4942 = t4561 * t4604;
t4567 = pkin(2) * cos(t4612) + t4630 * pkin(1);
t4562 = 0.1e1 / t4567;
t4941 = t4562 * t4622;
t4940 = t4562 * t4631;
t4563 = 0.1e1 / t4567 ^ 2;
t4597 = t4622 ^ 2;
t4939 = t4563 * t4597;
t4607 = t4631 ^ 2;
t4938 = t4563 * t4607;
t4568 = pkin(2) * cos(t4613) + t4633 * pkin(1);
t4564 = 0.1e1 / t4568;
t4937 = t4564 * t4625;
t4936 = t4564 * t4634;
t4565 = 0.1e1 / t4568 ^ 2;
t4601 = t4625 ^ 2;
t4935 = t4565 * t4601;
t4610 = t4634 ^ 2;
t4934 = t4565 * t4610;
t4602 = t4626 ^ 2;
t4569 = pkin(1) * t4626 + t4602 * t4990 - pkin(2);
t4933 = t4569 * t4618;
t4932 = t4569 * t4619;
t4605 = t4629 ^ 2;
t4570 = pkin(1) * t4629 + t4605 * t4990 - pkin(2);
t4931 = t4570 * t4621;
t4930 = t4570 * t4622;
t4608 = t4632 ^ 2;
t4571 = pkin(1) * t4632 + t4608 * t4990 - pkin(2);
t4929 = t4571 * t4624;
t4928 = t4571 * t4625;
t4927 = (pkin(1) + 0.2e1 * t4985) * t4617;
t4926 = (pkin(1) + 0.2e1 * t4984) * t4620;
t4925 = (pkin(1) + 0.2e1 * t4983) * t4623;
t4924 = t4584 * t4628;
t4923 = t4585 * t4631;
t4922 = t4586 * t4634;
t4921 = t4587 * t4628;
t4920 = t4588 * t4631;
t4919 = t4589 * t4634;
t4918 = t4590 * t4619;
t4917 = t4594 * t4622;
t4916 = t4598 * t4625;
t4913 = t4618 * t4619;
t4912 = t4618 * t4627;
t4909 = t4621 * t4622;
t4908 = t4621 * t4630;
t4905 = t4624 * t4625;
t4904 = t4624 * t4633;
t4903 = t4626 * t4617;
t4902 = t4626 * t4627;
t4900 = t4629 * t4620;
t4899 = t4629 * t4630;
t4897 = t4632 * t4623;
t4896 = t4632 * t4633;
t4645 = 1 / pkin(2);
t4894 = t4645 * t4646;
t4893 = 2 * t4973;
t4892 = 2 * t4972;
t4891 = -0.2e1 * t4914;
t4890 = -0.2e1 * t4910;
t4889 = -0.2e1 * t4906;
t4888 = -0.2e1 * t4902;
t4887 = -0.2e1 * t4899;
t4886 = -0.2e1 * t4896;
t4885 = pkin(2) * t4903;
t4884 = pkin(2) * t4900;
t4883 = pkin(2) * t4897;
t4882 = pkin(4) * t4560 * t4627;
t4881 = pkin(4) * t4562 * t4630;
t4880 = pkin(4) * t4564 * t4633;
t4527 = -t4901 * t4915 + (t4602 - 0.1e1) * t4619 * pkin(2);
t4603 = t4627 ^ 2;
t4822 = t4617 * t4913;
t4671 = pkin(1) * t4822 + (t4822 * t4990 + t4901) * t4626;
t4470 = (-t4584 * t4932 + t4587 * t4927) * t4603 + (t4584 * t4671 + t4587 * t4933) * t4627 + t4527 * t4584 - t4587 * t4885;
t4879 = t4470 * t4951;
t4471 = (t4584 * t4927 + t4587 * t4932) * t4603 + (t4584 * t4933 - t4587 * t4671) * t4627 - t4527 * t4587 - t4584 * t4885;
t4878 = t4471 * t4951;
t4528 = -t4898 * t4911 + (t4605 - 0.1e1) * t4622 * pkin(2);
t4606 = t4630 ^ 2;
t4821 = t4620 * t4909;
t4670 = pkin(1) * t4821 + (t4821 * t4990 + t4898) * t4629;
t4472 = (-t4585 * t4930 + t4588 * t4926) * t4606 + (t4585 * t4670 + t4588 * t4931) * t4630 + t4528 * t4585 - t4588 * t4884;
t4877 = t4472 * t4949;
t4473 = (t4585 * t4926 + t4588 * t4930) * t4606 + (t4585 * t4931 - t4588 * t4670) * t4630 - t4528 * t4588 - t4585 * t4884;
t4876 = t4473 * t4949;
t4529 = -t4895 * t4907 + (t4608 - 0.1e1) * t4625 * pkin(2);
t4609 = t4633 ^ 2;
t4820 = t4623 * t4905;
t4669 = pkin(1) * t4820 + (t4820 * t4990 + t4895) * t4632;
t4474 = (-t4586 * t4928 + t4589 * t4925) * t4609 + (t4586 * t4669 + t4589 * t4929) * t4633 + t4529 * t4586 - t4589 * t4883;
t4875 = t4474 * t4947;
t4475 = (t4586 * t4925 + t4589 * t4928) * t4609 + (t4586 * t4929 - t4589 * t4669) * t4633 - t4529 * t4589 - t4586 * t4883;
t4874 = t4475 * t4947;
t4819 = t4963 / 0.2e1;
t4485 = t4646 * t4819;
t4825 = t4590 * t4894;
t4476 = t4512 * t4825 + t4485;
t4873 = t4476 * t4951;
t4818 = t4962 / 0.2e1;
t4486 = t4646 * t4818;
t4824 = t4594 * t4894;
t4477 = t4513 * t4824 + t4486;
t4872 = t4477 * t4949;
t4817 = t4961 / 0.2e1;
t4487 = t4646 * t4817;
t4823 = t4598 * t4894;
t4478 = t4514 * t4823 + t4487;
t4871 = t4478 * t4947;
t4870 = t4536 * t4918;
t4869 = t4646 * t4951;
t4868 = t4538 * t4917;
t4867 = t4646 * t4949;
t4866 = t4540 * t4916;
t4865 = t4646 * t4947;
t4554 = t4902 - t4915;
t4864 = t4554 * t4945;
t4555 = t4899 - t4911;
t4863 = t4555 * t4941;
t4556 = t4896 - t4907;
t4862 = t4556 * t4937;
t4557 = t4618 * t4626 + t4914;
t4861 = t4557 * t4945;
t4558 = t4621 * t4629 + t4910;
t4860 = t4558 * t4941;
t4559 = t4624 * t4632 + t4906;
t4859 = t4559 * t4937;
t4858 = t4560 * t4924;
t4857 = t4560 * t4921;
t4856 = t4560 * t4918;
t4855 = t4590 * t4944;
t4854 = t4560 * t4913;
t4853 = t4618 * t4944;
t4578 = t4584 ^ 2;
t4852 = t4578 * t4942;
t4581 = t4587 ^ 2;
t4851 = t4581 * t4942;
t4592 = t4618 ^ 2;
t4850 = t4592 * t4942;
t4849 = t4561 * t4912;
t4848 = t4561 * t4619 * t4628;
t4847 = t4562 * t4923;
t4846 = t4562 * t4920;
t4845 = t4562 * t4917;
t4844 = t4594 * t4940;
t4843 = t4562 * t4909;
t4842 = t4621 * t4940;
t4579 = t4585 ^ 2;
t4841 = t4579 * t4938;
t4582 = t4588 ^ 2;
t4840 = t4582 * t4938;
t4596 = t4621 ^ 2;
t4839 = t4596 * t4938;
t4838 = t4563 * t4908;
t4837 = t4563 * t4622 * t4631;
t4836 = t4564 * t4922;
t4835 = t4564 * t4919;
t4834 = t4564 * t4916;
t4833 = t4598 * t4936;
t4832 = t4564 * t4905;
t4831 = t4624 * t4936;
t4580 = t4586 ^ 2;
t4830 = t4580 * t4934;
t4583 = t4589 ^ 2;
t4829 = t4583 * t4934;
t4600 = t4624 ^ 2;
t4828 = t4600 * t4934;
t4827 = t4565 * t4904;
t4826 = t4565 * t4625 * t4634;
t4816 = t4894 / 0.2e1;
t4815 = t4560 * t4603 * t4991;
t4814 = t4562 * t4606 * t4991;
t4813 = t4564 * t4609 * t4991;
t4812 = t4619 * t4882;
t4811 = t4628 * t4882;
t4810 = t4622 * t4881;
t4809 = t4631 * t4881;
t4808 = t4625 * t4880;
t4807 = t4634 * t4880;
t4806 = pkin(4) * t4854;
t4805 = pkin(4) * t4843;
t4804 = pkin(4) * t4832;
t4803 = t4945 * t4963;
t4802 = t4941 * t4962;
t4801 = t4937 * t4961;
t4800 = t4554 * t4858;
t4799 = t4554 * t4857;
t4798 = t4554 * t4856;
t4797 = t4555 * t4847;
t4796 = t4555 * t4846;
t4795 = t4555 * t4845;
t4794 = t4556 * t4836;
t4793 = t4556 * t4835;
t4792 = t4556 * t4834;
t4791 = t4557 * t4858;
t4790 = t4557 * t4857;
t4789 = t4557 * t4856;
t4788 = t4558 * t4847;
t4787 = t4558 * t4846;
t4786 = t4558 * t4845;
t4785 = t4559 * t4836;
t4784 = t4559 * t4835;
t4783 = t4559 * t4834;
t4782 = t4584 * t4855;
t4781 = t4584 * t4853;
t4780 = t4587 * t4855;
t4779 = t4587 * t4853;
t4778 = t4942 * t4995;
t4777 = t4584 * t4848;
t4776 = t4587 * t4848;
t4775 = t4592 * t4848;
t4774 = t4604 * t4849;
t4773 = t4585 * t4844;
t4772 = t4585 * t4842;
t4771 = t4588 * t4844;
t4770 = t4588 * t4842;
t4769 = t4938 * t4994;
t4768 = t4585 * t4837;
t4767 = t4588 * t4837;
t4766 = t4596 * t4837;
t4765 = t4607 * t4838;
t4764 = t4586 * t4833;
t4763 = t4586 * t4831;
t4762 = t4589 * t4833;
t4761 = t4589 * t4831;
t4760 = t4934 * t4993;
t4759 = t4586 * t4826;
t4758 = t4589 * t4826;
t4757 = t4600 * t4826;
t4756 = t4610 * t4827;
t4755 = t4590 * t4816;
t4754 = t4594 * t4816;
t4753 = t4598 * t4816;
t4752 = t4628 * t4815;
t4751 = t4631 * t4814;
t4750 = t4634 * t4813;
t4749 = t4617 * t4811;
t4748 = t4626 * t4811;
t4747 = t4620 * t4809;
t4746 = t4629 * t4809;
t4745 = t4623 * t4807;
t4744 = t4632 * t4807;
t4743 = pkin(4) * t4781;
t4742 = pkin(4) * t4779;
t4741 = pkin(4) * t4772;
t4740 = pkin(4) * t4770;
t4739 = pkin(4) * t4763;
t4738 = pkin(4) * t4761;
t4737 = t4536 * t4798;
t4736 = t4536 * t4789;
t4735 = t4538 * t4795;
t4734 = t4538 * t4786;
t4733 = t4540 * t4792;
t4732 = t4540 * t4783;
t4731 = t4554 * t4782;
t4730 = t4554 * t4780;
t4729 = t4555 * t4773;
t4728 = t4555 * t4771;
t4727 = t4556 * t4764;
t4726 = t4556 * t4762;
t4725 = t4557 * t4782;
t4724 = t4557 * t4780;
t4723 = t4558 * t4773;
t4722 = t4558 * t4771;
t4721 = t4559 * t4764;
t4720 = t4559 * t4762;
t4719 = t4848 * t4912;
t4718 = t4837 * t4908;
t4717 = t4826 * t4904;
t4716 = t4819 * t4951;
t4715 = -t4803 / 0.2e1;
t4714 = t4818 * t4949;
t4713 = -t4802 / 0.2e1;
t4712 = t4817 * t4947;
t4711 = -t4801 / 0.2e1;
t4710 = -t4924 * t4963 / 0.2e1;
t4709 = -t4923 * t4962 / 0.2e1;
t4708 = -t4922 * t4961 / 0.2e1;
t4707 = t4819 * t4921;
t4706 = t4818 * t4920;
t4705 = t4817 * t4919;
t4704 = t4587 * t4749;
t4703 = t4587 * t4748;
t4702 = t4588 * t4747;
t4701 = t4588 * t4746;
t4700 = t4589 * t4745;
t4699 = t4589 * t4744;
t4695 = t4470 * t4536 * t4782;
t4694 = t4471 * t4536 * t4780;
t4693 = t4472 * t4538 * t4773;
t4692 = t4473 * t4538 * t4771;
t4691 = t4474 * t4540 * t4764;
t4690 = t4475 * t4540 * t4762;
t4689 = t4536 * t4731;
t4688 = t4536 * t4730;
t4687 = t4536 * t4725;
t4686 = t4536 * t4724;
t4685 = t4538 * t4729;
t4684 = t4538 * t4728;
t4683 = t4538 * t4723;
t4682 = t4538 * t4722;
t4681 = t4540 * t4727;
t4680 = t4540 * t4726;
t4679 = t4540 * t4721;
t4678 = t4540 * t4720;
t4677 = t4560 * t4710;
t4676 = t4560 * t4707;
t4675 = t4562 * t4709;
t4674 = t4562 * t4706;
t4673 = t4564 * t4708;
t4672 = t4564 * t4705;
t4668 = t4476 * t4982 + t4619 * t4815;
t4667 = t4477 * t4981 + t4622 * t4814;
t4666 = t4478 * t4980 + t4625 * t4813;
t4509 = (t4602 - 0.1e1 / 0.2e1) * t4912 + (t4603 - 0.1e1 / 0.2e1) * t4903;
t4510 = (t4605 - 0.1e1 / 0.2e1) * t4908 + (t4606 - 0.1e1 / 0.2e1) * t4900;
t4511 = (t4608 - 0.1e1 / 0.2e1) * t4904 + (t4609 - 0.1e1 / 0.2e1) * t4897;
t4551 = t4557 ^ 2;
t4552 = t4558 ^ 2;
t4553 = t4559 ^ 2;
t4654 = t4540 * (t4474 * t4589 - t4475 * t4586) * t4833;
t4655 = t4538 * (t4472 * t4588 - t4473 * t4585) * t4844;
t4656 = t4536 * (t4470 * t4587 - t4471 * t4584) * t4855;
t4665 = (-t4618 * t4656 - t4621 * t4655 - t4624 * t4654) * t4973 + (-t4627 * t4656 - t4630 * t4655 - t4633 * t4654) * t4972 + (t4470 * t4471 * t4950 + t4472 * t4473 * t4948 + t4474 * t4475 * t4946) * t4971 + (-t4509 * t4778 - t4510 * t4769 - t4511 * t4760) * t4988 + (-t4551 * t4778 - t4552 * t4769 - t4553 * t4760) * MDP(11) + (-t4756 * t4993 - t4765 * t4994 - t4774 * t4995) * t4989 + (-t4592 * t4778 - t4596 * t4769 - t4600 * t4760) * MDP(4) + (-t4760 - t4769 - t4778) * MDP(1);
t4649 = t4564 * (t4474 * t4866 + t4708);
t4651 = t4562 * (t4472 * t4868 + t4709);
t4653 = t4560 * (t4470 * t4870 + t4710);
t4664 = (t4618 * t4653 + t4621 * t4651 + t4624 * t4649) * t4973 + (t4627 * t4653 + t4630 * t4651 + t4633 * t4649) * t4972 + (-t4470 * t4716 - t4472 * t4714 - t4474 * t4712) * t4971 + (t4509 * t4777 + t4510 * t4768 + t4511 * t4759) * t4988 + (t4551 * t4777 + t4552 * t4768 + t4553 * t4759) * MDP(11) + (t4584 * t4719 + t4585 * t4718 + t4586 * t4717) * t4989 + (t4584 * t4775 + t4585 * t4766 + t4586 * t4757) * MDP(4) + (t4759 + t4768 + t4777) * MDP(1);
t4648 = t4564 * (t4475 * t4866 + t4705);
t4650 = t4562 * (t4473 * t4868 + t4706);
t4652 = t4560 * (t4471 * t4870 + t4707);
t4663 = (t4618 * t4652 + t4621 * t4650 + t4624 * t4648) * t4973 + (t4627 * t4652 + t4630 * t4650 + t4633 * t4648) * t4972 + (-t4471 * t4716 - t4473 * t4714 - t4475 * t4712) * t4971 + (-t4509 * t4776 - t4510 * t4767 - t4511 * t4758) * t4988 + (-t4551 * t4776 - t4552 * t4767 - t4553 * t4758) * MDP(11) + (-t4587 * t4719 - t4588 * t4718 - t4589 * t4717) * t4989 + (-t4587 * t4775 - t4588 * t4766 - t4589 * t4757) * MDP(4) + (-t4758 - t4767 - t4776) * MDP(1);
t4450 = t4470 * t4869;
t4419 = t4500 * t4825 - t4450;
t4662 = t4419 * t4982 + t4584 * t4752;
t4451 = t4471 * t4869;
t4420 = t4501 * t4825 - t4451;
t4661 = -t4420 * t4982 + t4587 * t4752;
t4452 = t4472 * t4867;
t4421 = t4502 * t4824 - t4452;
t4660 = t4421 * t4981 + t4585 * t4751;
t4453 = t4473 * t4867;
t4422 = t4503 * t4824 - t4453;
t4659 = -t4422 * t4981 + t4588 * t4751;
t4454 = t4474 * t4865;
t4423 = t4504 * t4823 - t4454;
t4658 = t4423 * t4980 + t4586 * t4750;
t4455 = t4475 * t4865;
t4424 = t4505 * t4823 - t4455;
t4657 = -t4424 * t4980 + t4589 * t4750;
t4535 = t4632 * t4808;
t4534 = t4629 * t4810;
t4533 = t4626 * t4812;
t4532 = t4623 * t4808;
t4531 = t4620 * t4810;
t4530 = t4617 * t4812;
t4520 = t4586 * t4744;
t4519 = t4585 * t4746;
t4518 = t4584 * t4748;
t4517 = t4586 * t4745;
t4516 = t4585 * t4747;
t4515 = t4584 * t4749;
t4481 = t4804 + t4817;
t4480 = t4805 + t4818;
t4479 = t4806 + t4819;
t4469 = -t4481 * t4623 + t4535;
t4468 = -t4480 * t4620 + t4534;
t4467 = -t4479 * t4617 + t4533;
t4466 = t4481 * t4632 + t4532;
t4465 = t4480 * t4629 + t4531;
t4464 = t4479 * t4626 + t4530;
t4463 = -pkin(1) * t4832 + t4478 * t4986;
t4462 = t4804 + (t4514 * t4753 + t4487) * t4991;
t4461 = -pkin(1) * t4843 + t4477 * t4986;
t4460 = t4805 + (t4513 * t4754 + t4486) * t4991;
t4459 = -pkin(1) * t4854 + t4476 * t4986;
t4458 = t4806 + (t4512 * t4755 + t4485) * t4991;
t4448 = -t4462 * t4623 + t4535;
t4447 = t4462 * t4632 + t4532;
t4446 = -t4460 * t4620 + t4534;
t4445 = t4460 * t4629 + t4531;
t4444 = -t4458 * t4617 + t4533;
t4443 = t4458 * t4626 + t4530;
t4442 = -t4738 - t4874;
t4441 = t4739 - t4875;
t4440 = -t4740 - t4876;
t4439 = t4741 - t4877;
t4438 = -t4742 - t4878;
t4437 = t4743 - t4879;
t4436 = -t4442 * t4623 - t4699;
t4435 = -t4441 * t4623 + t4520;
t4434 = -t4440 * t4620 - t4701;
t4433 = -t4439 * t4620 + t4519;
t4432 = -t4438 * t4617 - t4703;
t4431 = -t4437 * t4617 + t4518;
t4430 = t4442 * t4632 - t4700;
t4429 = t4441 * t4632 + t4517;
t4428 = t4440 * t4629 - t4702;
t4427 = t4439 * t4629 + t4516;
t4426 = t4438 * t4626 - t4704;
t4425 = t4437 * t4626 + t4515;
t4418 = pkin(1) * t4761 + t4424 * t4986;
t4417 = -pkin(1) * t4763 + t4423 * t4986;
t4416 = t4738 + (t4505 * t4753 - t4455) * t4992;
t4415 = t4739 + (t4504 * t4753 - t4454) * t4991;
t4414 = pkin(1) * t4770 + t4422 * t4986;
t4413 = -pkin(1) * t4772 + t4421 * t4986;
t4412 = t4740 + (t4503 * t4754 - t4453) * t4992;
t4411 = t4741 + (t4502 * t4754 - t4452) * t4991;
t4410 = pkin(1) * t4779 + t4420 * t4986;
t4409 = -pkin(1) * t4781 + t4419 * t4986;
t4408 = t4742 + (t4501 * t4755 - t4451) * t4992;
t4407 = t4743 + (t4500 * t4755 - t4450) * t4991;
t4406 = -t4416 * t4632 - t4700;
t4405 = t4416 * t4623 - t4699;
t4404 = -t4415 * t4623 + t4520;
t4403 = t4415 * t4632 + t4517;
t4402 = -t4412 * t4629 - t4702;
t4401 = t4412 * t4620 - t4701;
t4400 = -t4411 * t4620 + t4519;
t4399 = t4411 * t4629 + t4516;
t4398 = -t4408 * t4626 - t4704;
t4397 = t4408 * t4617 - t4703;
t4396 = -t4407 * t4617 + t4518;
t4395 = t4407 * t4626 + t4515;
t4394 = t4463 * t4886 + t4623 * t4666;
t4393 = t4463 * t4889 - t4632 * t4666;
t4392 = t4461 * t4887 + t4620 * t4667;
t4391 = t4461 * t4890 - t4629 * t4667;
t4390 = t4459 * t4888 + t4617 * t4668;
t4389 = t4459 * t4891 - t4626 * t4668;
t4388 = t4418 * t4886 - t4623 * t4657;
t4387 = t4417 * t4886 + t4623 * t4658;
t4386 = t4418 * t4889 + t4632 * t4657;
t4385 = t4417 * t4889 - t4632 * t4658;
t4384 = t4414 * t4887 - t4620 * t4659;
t4383 = t4413 * t4887 + t4620 * t4660;
t4382 = t4414 * t4890 + t4629 * t4659;
t4381 = t4413 * t4890 - t4629 * t4660;
t4380 = t4410 * t4888 - t4617 * t4661;
t4379 = t4409 * t4888 + t4617 * t4662;
t4378 = t4410 * t4891 + t4626 * t4661;
t4377 = t4409 * t4891 - t4626 * t4662;
t1 = [(t4829 + t4840 + t4851) * MDP(1) + (t4581 * t4850 + t4582 * t4839 + t4583 * t4828) * MDP(4) + (t4581 * t4774 + t4582 * t4765 + t4583 * t4756) * t4989 + (-t4618 * t4694 - t4621 * t4692 - t4624 * t4690) * t4893 + (-t4627 * t4694 - t4630 * t4692 - t4633 * t4690) * t4892 + (t4471 ^ 2 * t4950 + t4473 ^ 2 * t4948 + t4475 ^ 2 * t4946) * t4971 + (t4551 * t4851 + t4552 * t4840 + t4553 * t4829) * MDP(11) + (t4509 * t4851 + t4510 * t4840 + t4511 * t4829) * t4988 + (t4420 * t4790 + t4422 * t4787 + t4424 * t4784 + (-t4471 * t4686 - t4473 * t4682 - t4475 * t4678 + (t4501 * t4724 + t4503 * t4722 + t4505 * t4720) * t4645) * t4646) * MDP(13) + (t4420 * t4799 + t4422 * t4796 + t4424 * t4793 + (-t4471 * t4688 - t4473 * t4684 - t4475 * t4680 + (t4501 * t4730 + t4503 * t4728 + t4505 * t4726) * t4645) * t4646) * MDP(14) + (-t4420 * t4878 - t4422 * t4876 - t4424 * t4874 + (t4420 * t4959 + t4422 * t4957 + t4424 * t4955) * t4645) * t4964 + (t4378 * t4857 + t4382 * t4846 + t4386 * t4835 + (-t4398 * t4878 - t4402 * t4876 - t4406 * t4874 + (t4426 * t4959 + t4428 * t4957 + t4430 * t4955) * t4645) * t4646) * MDP(16) + (t4380 * t4857 + t4384 * t4846 + t4388 * t4835 + (-t4397 * t4878 - t4401 * t4876 - t4405 * t4874 + (t4432 * t4959 + t4434 * t4957 + t4436 * t4955) * t4645) * t4646) * MDP(17) + MDP(18); (t4419 * t4790 + t4421 * t4787 + t4423 * t4784 + (t4471 * t4687 + t4473 * t4683 + t4475 * t4679 + (-t4501 * t4725 - t4503 * t4723 - t4505 * t4721) * t4645) * t4646) * MDP(13) + (t4419 * t4799 + t4421 * t4796 + t4423 * t4793 + (t4471 * t4689 + t4473 * t4685 + t4475 * t4681 + (-t4501 * t4731 - t4503 * t4729 - t4505 * t4727) * t4645) * t4646) * MDP(14) + (-t4419 * t4878 - t4421 * t4876 - t4423 * t4874 + (t4419 * t4959 + t4421 * t4957 + t4423 * t4955) * t4645) * t4964 + (t4377 * t4857 + t4381 * t4846 + t4385 * t4835 + (-t4395 * t4878 - t4399 * t4876 - t4403 * t4874 + (t4425 * t4959 + t4427 * t4957 + t4429 * t4955) * t4645) * t4646) * MDP(16) + (t4379 * t4857 + t4383 * t4846 + t4387 * t4835 + (-t4396 * t4878 - t4400 * t4876 - t4404 * t4874 + (t4431 * t4959 + t4433 * t4957 + t4435 * t4955) * t4645) * t4646) * MDP(17) + t4665; (t4476 * t4790 + t4477 * t4787 + t4478 * t4784 + (t4471 * t4736 + t4473 * t4734 + t4475 * t4732 + (-t4501 * t4789 - t4503 * t4786 - t4505 * t4783) * t4645) * t4646) * MDP(13) + (t4476 * t4799 + t4477 * t4796 + t4478 * t4793 + (t4471 * t4737 + t4473 * t4735 + t4475 * t4733 + (-t4501 * t4798 - t4503 * t4795 - t4505 * t4792) * t4645) * t4646) * MDP(14) + (-t4471 * t4873 - t4473 * t4872 - t4475 * t4871 + (t4476 * t4959 + t4477 * t4957 + t4478 * t4955) * t4645) * t4964 + (t4389 * t4857 + t4391 * t4846 + t4393 * t4835 + (-t4443 * t4878 - t4445 * t4876 - t4447 * t4874 + (t4464 * t4959 + t4465 * t4957 + t4466 * t4955) * t4645) * t4646) * MDP(16) + (t4390 * t4857 + t4392 * t4846 + t4394 * t4835 + (-t4444 * t4878 - t4446 * t4876 - t4448 * t4874 + (t4467 * t4959 + t4468 * t4957 + t4469 * t4955) * t4645) * t4646) * MDP(17) + t4663; (-t4420 * t4791 - t4422 * t4788 - t4424 * t4785 + (-t4470 * t4686 - t4472 * t4682 - t4474 * t4678 + (t4500 * t4724 + t4502 * t4722 + t4504 * t4720) * t4645) * t4646) * MDP(13) + (-t4420 * t4800 - t4422 * t4797 - t4424 * t4794 + (-t4470 * t4688 - t4472 * t4684 - t4474 * t4680 + (t4500 * t4730 + t4502 * t4728 + t4504 * t4726) * t4645) * t4646) * MDP(14) + (-t4420 * t4879 - t4422 * t4877 - t4424 * t4875 + (t4420 * t4960 + t4422 * t4958 + t4424 * t4956) * t4645) * t4964 + (-t4378 * t4858 - t4382 * t4847 - t4386 * t4836 + (-t4398 * t4879 - t4402 * t4877 - t4406 * t4875 + (t4426 * t4960 + t4428 * t4958 + t4430 * t4956) * t4645) * t4646) * MDP(16) + (-t4380 * t4858 - t4384 * t4847 - t4388 * t4836 + (-t4397 * t4879 - t4401 * t4877 - t4405 * t4875 + (t4432 * t4960 + t4434 * t4958 + t4436 * t4956) * t4645) * t4646) * MDP(17) + t4665; (t4830 + t4841 + t4852) * MDP(1) + (t4578 * t4850 + t4579 * t4839 + t4580 * t4828) * MDP(4) + (t4578 * t4774 + t4579 * t4765 + t4580 * t4756) * t4989 + (t4618 * t4695 + t4621 * t4693 + t4624 * t4691) * t4893 + (t4627 * t4695 + t4630 * t4693 + t4633 * t4691) * t4892 + (t4470 ^ 2 * t4950 + t4472 ^ 2 * t4948 + t4474 ^ 2 * t4946) * t4971 + (t4551 * t4852 + t4552 * t4841 + t4553 * t4830) * MDP(11) + (t4509 * t4852 + t4510 * t4841 + t4511 * t4830) * t4988 + (-t4419 * t4791 - t4421 * t4788 - t4423 * t4785 + (t4470 * t4687 + t4472 * t4683 + t4474 * t4679 + (-t4500 * t4725 - t4502 * t4723 - t4504 * t4721) * t4645) * t4646) * MDP(13) + (-t4419 * t4800 - t4421 * t4797 - t4423 * t4794 + (t4470 * t4689 + t4472 * t4685 + t4474 * t4681 + (-t4500 * t4731 - t4502 * t4729 - t4504 * t4727) * t4645) * t4646) * MDP(14) + (-t4419 * t4879 - t4421 * t4877 - t4423 * t4875 + (t4419 * t4960 + t4421 * t4958 + t4423 * t4956) * t4645) * t4964 + (-t4377 * t4858 - t4381 * t4847 - t4385 * t4836 + (-t4395 * t4879 - t4399 * t4877 - t4403 * t4875 + (t4425 * t4960 + t4427 * t4958 + t4429 * t4956) * t4645) * t4646) * MDP(16) + (-t4379 * t4858 - t4383 * t4847 - t4387 * t4836 + (-t4396 * t4879 - t4400 * t4877 - t4404 * t4875 + (t4431 * t4960 + t4433 * t4958 + t4435 * t4956) * t4645) * t4646) * MDP(17) + MDP(18); (-t4476 * t4791 - t4477 * t4788 - t4478 * t4785 + (t4470 * t4736 + t4472 * t4734 + t4474 * t4732 + (-t4500 * t4789 - t4502 * t4786 - t4504 * t4783) * t4645) * t4646) * MDP(13) + (-t4476 * t4800 - t4477 * t4797 - t4478 * t4794 + (t4470 * t4737 + t4472 * t4735 + t4474 * t4733 + (-t4500 * t4798 - t4502 * t4795 - t4504 * t4792) * t4645) * t4646) * MDP(14) + (-t4470 * t4873 - t4472 * t4872 - t4474 * t4871 + (t4476 * t4960 + t4477 * t4958 + t4478 * t4956) * t4645) * t4964 + (-t4389 * t4858 - t4391 * t4847 - t4393 * t4836 + (-t4443 * t4879 - t4445 * t4877 - t4447 * t4875 + (t4464 * t4960 + t4465 * t4958 + t4466 * t4956) * t4645) * t4646) * MDP(16) + (-t4390 * t4858 - t4392 * t4847 - t4394 * t4836 + (-t4444 * t4879 - t4446 * t4877 - t4448 * t4875 + (t4467 * t4960 + t4468 * t4958 + t4469 * t4956) * t4645) * t4646) * MDP(17) + t4664; (-t4420 * t4861 - t4422 * t4860 - t4424 * t4859 + (t4559 * t4672 + t4558 * t4674 + t4557 * t4676 + (t4512 * t4724 + t4513 * t4722 + t4514 * t4720) * t4645) * t4646) * MDP(13) + (-t4420 * t4864 - t4422 * t4863 - t4424 * t4862 + (t4556 * t4672 + t4555 * t4674 + t4554 * t4676 + (t4512 * t4730 + t4513 * t4728 + t4514 * t4726) * t4645) * t4646) * MDP(14) + (t4424 * t4817 + t4422 * t4818 + t4420 * t4819 + (t4420 * t4954 + t4422 * t4953 + t4424 * t4952) * t4645) * t4964 + (-t4378 * t4945 - t4382 * t4941 - t4386 * t4937 + (t4406 * t4817 + t4402 * t4818 + t4398 * t4819 + (t4426 * t4954 + t4428 * t4953 + t4430 * t4952) * t4645) * t4646) * MDP(16) + (-t4380 * t4945 - t4384 * t4941 - t4388 * t4937 + (t4405 * t4817 + t4401 * t4818 + t4397 * t4819 + (t4432 * t4954 + t4434 * t4953 + t4436 * t4952) * t4645) * t4646) * MDP(17) + t4663; (-t4419 * t4861 - t4421 * t4860 - t4423 * t4859 + (t4559 * t4673 + t4558 * t4675 + t4557 * t4677 + (-t4512 * t4725 - t4513 * t4723 - t4514 * t4721) * t4645) * t4646) * MDP(13) + (-t4419 * t4864 - t4421 * t4863 - t4423 * t4862 + (t4556 * t4673 + t4555 * t4675 + t4554 * t4677 + (-t4512 * t4731 - t4513 * t4729 - t4514 * t4727) * t4645) * t4646) * MDP(14) + (t4423 * t4817 + t4421 * t4818 + t4419 * t4819 + (t4419 * t4954 + t4421 * t4953 + t4423 * t4952) * t4645) * t4964 + (-t4377 * t4945 - t4381 * t4941 - t4385 * t4937 + (t4403 * t4817 + t4399 * t4818 + t4395 * t4819 + (t4425 * t4954 + t4427 * t4953 + t4429 * t4952) * t4645) * t4646) * MDP(16) + (-t4379 * t4945 - t4383 * t4941 - t4387 * t4937 + (t4404 * t4817 + t4400 * t4818 + t4396 * t4819 + (t4431 * t4954 + t4433 * t4953 + t4435 * t4952) * t4645) * t4646) * MDP(17) + t4664; (t4935 + t4939 + t4943) * MDP(1) + (t4592 * t4943 + t4596 * t4939 + t4600 * t4935) * MDP(4) + (t4593 * t4849 + t4597 * t4838 + t4601 * t4827) * t4989 + (-t4618 * t4803 - t4621 * t4802 - t4624 * t4801) * t4973 + (-t4627 * t4803 - t4630 * t4802 - t4633 * t4801) * t4972 + (t4497 ^ 2 / t4526 ^ 2 / 0.4e1 + t4496 ^ 2 / t4525 ^ 2 / 0.4e1 + t4495 ^ 2 / t4524 ^ 2 / 0.4e1) * t4971 + (t4551 * t4943 + t4552 * t4939 + t4553 * t4935) * MDP(11) + (t4509 * t4943 + t4510 * t4939 + t4511 * t4935) * t4988 + (-t4476 * t4861 - t4477 * t4860 - t4478 * t4859 + (t4559 * t4711 + t4558 * t4713 + t4557 * t4715 + (-t4512 * t4789 - t4513 * t4786 - t4514 * t4783) * t4645) * t4646) * MDP(13) + (-t4476 * t4864 - t4477 * t4863 - t4478 * t4862 + (t4556 * t4711 + t4555 * t4713 + t4554 * t4715 + (-t4512 * t4798 - t4513 * t4795 - t4514 * t4792) * t4645) * t4646) * MDP(14) + (t4478 * t4817 + t4477 * t4818 + t4476 * t4819 + (t4476 * t4954 + t4477 * t4953 + t4478 * t4952) * t4645) * t4964 + (-t4389 * t4945 - t4391 * t4941 - t4393 * t4937 + (t4447 * t4817 + t4445 * t4818 + t4443 * t4819 + (t4464 * t4954 + t4465 * t4953 + t4466 * t4952) * t4645) * t4646) * MDP(16) + (-t4390 * t4945 - t4392 * t4941 - t4394 * t4937 + (t4448 * t4817 + t4446 * t4818 + t4444 * t4819 + (t4467 * t4954 + t4468 * t4953 + t4469 * t4952) * t4645) * t4646) * MDP(17) + MDP(18);];
%% Postprocessing: Reshape Output
% From vec2mat_3_matlab.m
res = [t1(1), t1(2), t1(3); t1(4), t1(5), t1(6); t1(7), t1(8), t1(9);];
MMX  = res;
