# Việc cần làm

1. Cài R và các phiên bản gói nêu trong README; smoke test Mỹ/Trung Quốc/Canada với `WTC_CAL=25 WTC_TEST=49`, sau đó chạy lại pilot với 250/999; so khớp cửa sổ mẫu và chỉ số cụm/episode với báo cáo trước. Sự khác nhau về phương pháp định nghĩa cluster phải được ghi rõ, không ép khớp p.
2. Chạy mặc định 18 nước với null phương sai cục bộ, lưu `results/core/summary.csv`, checkpoint và `sessionInfo.txt`; ghi ngày, phiên bản dữ liệu và thời gian chạy. Báo kết quả mới sau khi có đầu ra, không chuyển số pilot sang 18 nước.
3. Chạy `WTC_NULL=ar1`, `WTC_NULL=fourier` vào hai thư mục riêng, so band/cluster/episode, BH/BY/Holm. Kiểm tra null có tái tạo ACF, phổ, biên độ phương sai thay đổi theo thời gian và các đợt đứt gãy; kiểm tra kích thước và power bằng mô phỏng thiết kế trước.
4. Sửa thiếu CPI gốc Mỹ 2025-10 bằng đúng quan sát BIS cùng đơn vị nếu có, không tự nội suy. Cập nhật mùa vụ STL hai phía và kiểm tra hồi cứu/rolling endpoint; phân tích thêm EPU có mùa vụ và độ nhạy biến đổi MoM.
5. Chạy UCT và Brent trên dữ liệu mới tải, ghi mapping cột chính xác và biến đổi UCT mức/difference; chẩn đoán trị riêng/điều kiện của ma trận coherency; nếu muốn kiểm định partial/multiple phải viết surrogate chung giữ cấu trúc giữa controls và hiệu chuẩn riêng.
6. Ghi phiên bản upstream, điều khoản dữ liệu và sai khác gần cuối mẫu. Thiết kế xác nhận ngoài mẫu hoặc thêm biến kiểm soát theo lý thuyết. Diễn giải biến cố cùng tài liệu lịch sử và phương pháp nhận dạng vĩ mô độc lập trước khi nói đến tác động chính sách.
