# Kế hoạch nghiên cứu ban đầu và thay đổi

1. Ghép EPU quốc gia và CPI tháng, giữ lịch thật; tính EPU log MoM và INF log MoM từ CPI đã điều chỉnh mùa vụ. So sánh các cách điều chỉnh trước khi xem WTC.
2. Chạy Morlet CWT/WTC cho 18 nước, báo cáo ba nhóm chu kỳ 2–<8, 8–<32, ≥32 tháng; loại vùng ngoài COI, kiểm tra pha, vùng liên tục và các thay đổi ước lượng.
3. Kiểm tra AR(1), Fourier, BH/BY/Holm, thư viện khác, mùa vụ, độ dừng, lưới và làm trơn. Tách quan sát đồ thị khỏi phát biểu thống kê và nhân quả.
4. Sau khi bản đầu nhạy mạnh với smoother và kiểm định từng ô không bền, thử riêng Mỹ, Trung Quốc, Canada: smoother 5/0,25 và 9/1; null AR(1), Fourier, AR(1) phương sai biến thiên; kiểm định mean nhóm chu kỳ, max cluster mass và max episode mean. Bổ sung UCT/dầu cho Mỹ–Trung.
5. Đóng băng đặc tả 9/1 làm đường chạy lõi cho 18 nước. Đưa null thay thế, mùa vụ, thư viện và đặc tả 5/0,25 vào robustness; để conditional WTC trong `extensions/uct/`.

Kế hoạch ở đây ghi lại trình tự nghiên cứu, không phải đăng ký trước độc lập. Việc chọn smoother 9/1 dựa trên pilot đã xem dữ liệu, nên không mô tả phép thử cuối như xác nhận ngoài mẫu. Holm không bị xóa vì có p thô gần 0,05; đây là phép đối chiếu bảo thủ cạnh BH/BY. Khi công bố giả thuyết phải nêu trước family, thống kê và mức lỗi.
